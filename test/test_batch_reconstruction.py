"""Independent reconstruction oracles and native-state regression checks."""
from fractions import Fraction
import unittest
import flint
from aim044_lift.reconstruct import _library, recover_euclid, recover_many_flint


class BatchReconstructionTests(unittest.TestCase):
    def test_mixed_signed_batch_matches_euclidean_reference(self):
        modulus, height = 1 << 512, 1 << 200
        expected = [Fraction(0), Fraction(-17, 9), Fraction(23, 31), Fraction(2**150+1, 2**170+3)]
        residues = [x.numerator * pow(x.denominator, -1, modulus) % modulus for x in expected]
        result = recover_many_flint(iter(residues), modulus, height)
        self.assertEqual(result, [recover_euclid(r, modulus, height) for r in residues])
        self.assertEqual([Fraction(*p) for p in result], expected)

    def test_failure_does_not_poison_later_native_calls(self):
        with self.assertRaises(ValueError):
            recover_many_flint([0, 17], 65536, 8)
        self.assertEqual(recover_many_flint([0, 1], 65536, 8), [(0, 1), (1, 1)])

    def test_backend_initialization_preserves_process_thread_setting(self):
        previous = flint.ctx.threads
        try:
            flint.ctx.threads = 2
            _library.cache_clear()
            self.assertEqual(recover_many_flint([0], 4096, 16), [(0, 1)])
            self.assertEqual(flint.ctx.threads, 2)
        finally:
            flint.ctx.threads = previous
