"""Bounded rational recovery with a pinned FLINT half-GCD backend.

The classical Euclidean implementation is an independent reference only.
It does not have the fast reconstruction complexity used by the theorem.
"""

import ctypes
import importlib.metadata
import os
from functools import lru_cache
from math import gcd


def _check(residue, modulus, height):
    if modulus < 2 or height < 1 or not 0 <= residue < modulus:
        raise ValueError("Invalid reconstruction parameters")
    if modulus & (modulus - 1) or height & (height - 1):
        raise ValueError("This interface requires power-of-two modulus and height")
    if modulus.bit_length() <= 2 * height.bit_length():
        raise ValueError("The modulus must be strictly greater than twice height squared")


def recover_euclid(residue, modulus, height):
    """Slow, transparent reference returning a reduced numerator/denominator."""
    _check(residue, modulus, height)
    old_r, r = modulus, residue
    old_t, t = 0, 1
    while r > height:
        quotient, remainder = divmod(old_r, r)
        old_r, r = r, remainder
        old_t, t = t, old_t - quotient * t
    if t == 0:
        raise ValueError("No bounded reconstruction")
    numerator, denominator = (r, t) if t > 0 else (-r, -t)
    if (denominator > height or gcd(numerator, denominator) != 1
            or gcd(denominator, modulus) != 1
            or (numerator - residue * denominator) % modulus):
        raise ValueError("No bounded reconstruction")
    return numerator, denominator


class _Fmpq(ctypes.Structure):
    # FLINT fmpz is a signed machine word; the supported wheel ABI uses LP64.
    _fields_ = [("num", ctypes.c_long), ("den", ctypes.c_long)]


@lru_cache(maxsize=1)
def _library():
    if os.name != "posix" or ctypes.sizeof(ctypes.c_long) != 8:
        raise RuntimeError("The FLINT binding supports 64-bit Linux/macOS (LP64)")
    distribution = importlib.metadata.distribution("python-flint")
    if distribution.version != "0.8.0":
        raise RuntimeError("Install the pinned dependency python-flint==0.8.0")
    paths = [distribution.locate_file(item) for item in distribution.files or []
             if "libflint" in item.name
             and (".so" in item.name or item.name.endswith(".dylib"))]
    if len(paths) != 1:
        raise RuntimeError("Expected one bundled FLINT library in the pinned wheel")
    lib = ctypes.CDLL(str(paths[0]))
    word_pointer = ctypes.POINTER(ctypes.c_long)
    lib.fmpz_set_str.argtypes = [word_pointer, ctypes.c_char_p, ctypes.c_int]
    lib.fmpz_set_str.restype = ctypes.c_int
    lib.fmpz_get_str.argtypes = [ctypes.c_void_p, ctypes.c_int, word_pointer]
    lib.fmpz_get_str.restype = ctypes.c_void_p
    lib.fmpz_clear.argtypes = [word_pointer]
    lib.fmpz_clear.restype = None
    lib.flint_free.argtypes = [ctypes.c_void_p]
    lib.flint_free.restype = None
    lib.fmpq_reconstruct_fmpz_2.argtypes = [ctypes.POINTER(_Fmpq)] + [word_pointer] * 4
    lib.fmpq_reconstruct_fmpz_2.restype = ctypes.c_int
    return lib


def ensure_backend(name):
    """Check an explicitly selected backend before the costly lifting phase."""
    if name == "flint":
        _library()


def recover_flint(residue, modulus, height):
    """Recover with FLINT's recursive half-GCD bounded reconstruction.

    Base-16 conversion keeps data transfer linear in the binary bit length.
    The result is already reduced; no Python long-integer gcd is introduced.
    """
    return recover_many_flint((residue,), modulus, height)[0]


def recover_many_flint(residues, modulus, height):
    """Recover a batch with shared bounds and deterministic native cleanup.

    Bounds are converted once per solve. No process-global FLINT thread or
    precision setting is changed. The numerical reconstruction is unchanged.
    """
    _check(0, modulus, height)
    lib = _library()
    values = [ctypes.c_long(0) for _ in range(3)]
    result = _Fmpq(0, 1)
    pointers = [ctypes.pointer(value) for value in values]
    num_ptr = ctypes.cast(ctypes.byref(result), ctypes.POINTER(ctypes.c_long))
    den_ptr = ctypes.cast(ctypes.byref(result, ctypes.sizeof(ctypes.c_long)),
                          ctypes.POINTER(ctypes.c_long))
    try:
        for pointer, value in zip(pointers[1:], (modulus, height)):
            if lib.fmpz_set_str(pointer, format(value, "x").encode("ascii"), 16):
                raise RuntimeError("FLINT rejected a hexadecimal integer")
        output = []
        for residue in residues:
            _check(residue, modulus, height)
            if lib.fmpz_set_str(pointers[0], format(residue, "x").encode("ascii"), 16):
                raise RuntimeError("FLINT rejected a hexadecimal integer")
            if not lib.fmpq_reconstruct_fmpz_2(ctypes.byref(result), pointers[0],
                                             pointers[1], pointers[2], pointers[2]):
                raise ValueError("No bounded reconstruction")
            pair = []
            for pointer in (num_ptr, den_ptr):
                raw = lib.fmpz_get_str(None, 16, pointer)
                if not raw:
                    raise MemoryError("FLINT integer conversion failed")
                try:
                    pair.append(int(ctypes.string_at(raw), 16))
                finally:
                    lib.flint_free(raw)
            numerator, denominator = pair
            if abs(numerator) > height or not 1 <= denominator <= height or not denominator & 1:
                raise RuntimeError("FLINT returned an invalid bounded reconstruction")
            output.append((numerator, denominator))
        return output
    finally:
        for pointer in pointers + [num_ptr, den_ptr]:
            lib.fmpz_clear(pointer)
