import unittest
from tinygrad import Tensor
from fields import FieldElement, Field, Polynomial, xgcd

class TestFieldElement(unittest.TestCase):
    def setUp(self):
        self.field = Field(97)  # Using a small prime for easier testing

    def test_addition(self):
        a = FieldElement(10, self.field)
        b = FieldElement(20, self.field)
        self.assertEqual((a + b).value, 30)

    def test_subtraction(self):
        a = FieldElement(30, self.field)
        b = FieldElement(20, self.field)
        self.assertEqual((a - b).value, 10)

    def test_multiplication(self):
        a = FieldElement(10, self.field)
        b = FieldElement(20, self.field)
        self.assertEqual((a * b).value, 6)  # (10 * 20) % 97 = 6

    def test_division(self):
        a = FieldElement(10, self.field)
        b = FieldElement(20, self.field)
        self.assertEqual((a / b).value, 49)  # (10 * 20^-1) % 97 = 5

    def test_exponentiation(self):
        a = FieldElement(2, self.field)
        self.assertEqual((a ^ 5).value, 32)

    def test_equality(self):
        a = FieldElement(10, self.field)
        b = FieldElement(10, self.field)
        c = FieldElement(20, self.field)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)

class TestField(unittest.TestCase):
    def setUp(self):
        self.field = Field(97)

    def test_zero_and_one(self):
        self.assertEqual(self.field.zero().value, 0)
        self.assertEqual(self.field.one().value, 1)

    def test_addition(self):
        a = FieldElement(10, self.field)
        b = FieldElement(20, self.field)
        self.assertEqual(self.field.add(a, b).value, 30)

    def test_subtraction(self):
        a = FieldElement(30, self.field)
        b = FieldElement(20, self.field)
        self.assertEqual(self.field.sub(a, b).value, 10)

    def test_multiplication(self):
        a = FieldElement(10, self.field)
        b = FieldElement(20, self.field)
        self.assertEqual(self.field.mul(a, b).value, 6)

    def test_division(self):
        a = FieldElement(10, self.field)
        b = FieldElement(20, self.field)
        self.assertEqual(self.field.divide(a, b).value, 49)

    def test_inverse(self):
        a = FieldElement(10, self.field)
        inv_a = self.field.inverse(a)
        self.assertEqual((a * inv_a).value, 1)

    def test_negate(self):
        a = FieldElement(10, self.field)
        self.assertEqual(self.field.negate(a).value, 87)  # -10 mod 97 = 87

class TestPolynomial(unittest.TestCase):
    def setUp(self):
        self.field = Field(97)
        self.x = Polynomial([FieldElement(0, self.field), FieldElement(1, self.field)], self.field)  # x

    def test_addition(self):
        p1 = Polynomial([FieldElement(1, self.field), FieldElement(2, self.field)], self.field)  # 1 + 2x
        p2 = Polynomial([FieldElement(3, self.field), FieldElement(4, self.field)], self.field)  # 3 + 4x
        result = p1 + p2
        self.assertEqual([coeff.value for coeff in result.poly.tolist()], [4, 6])

    def test_subtraction(self):
        p1 = Polynomial([FieldElement(3, self.field), FieldElement(4, self.field)], self.field)  # 3 + 4x
        p2 = Polynomial([FieldElement(1, self.field), FieldElement(2, self.field)], self.field)  # 1 + 2x
        result = p1 - p2
        self.assertEqual([coeff.value for coeff in result.poly.tolist()], [2, 2])

    def test_multiplication(self):
        p1 = Polynomial([FieldElement(1, self.field), FieldElement(2, self.field)], self.field)  # 1 + 2x
        p2 = Polynomial([FieldElement(3, self.field), FieldElement(4, self.field)], self.field)  # 3 + 4x
        result = p1 * p2
        self.assertEqual([coeff.value for coeff in result.poly.tolist()], [3, 10, 8])

    def test_division(self):
        p1 = Polynomial([FieldElement(1, self.field), FieldElement(2, self.field), FieldElement(1, self.field)], self.field)  # 1 + 2x + x^2
        p2 = Polynomial([FieldElement(1, self.field), FieldElement(1, self.field)], self.field)  # 1 + x
        quotient, remainder = Polynomial.divide(p1, p2, self.field)
        self.assertEqual([coeff.value for coeff in quotient.poly.tolist()], [1, 1])
        self.assertEqual([coeff.value for coeff in remainder.poly.tolist()], [0])

    def test_evaluation(self):
        p = Polynomial([FieldElement(1, self.field), FieldElement(2, self.field), FieldElement(1, self.field)], self.field)  # 1 + 2x + x^2
        result = p.evaluate(FieldElement(2, self.field))
        self.assertEqual(result.value, 9)  # 1 + 2*2 + 2^2 = 9

    def test_interpolation(self):
        domain = [FieldElement(0, self.field), FieldElement(1, self.field), FieldElement(2, self.field)]
        values = [FieldElement(1, self.field), FieldElement(2, self.field), FieldElement(5, self.field)]
        p = Polynomial.interpolate_domain(domain, values, self.field)
        self.assertEqual([coeff.value for coeff in p.poly.tolist()], [1, 0, 1])  # 1 + x^2

    def test_zeroifier_domain(self):
        domain = [FieldElement(0, self.field), FieldElement(1, self.field), FieldElement(2, self.field)]
        z = Polynomial.zeroifier_domain(domain, self.field)
        self.assertEqual([coeff.value for coeff in z.poly.tolist()], [0, 3, 94, 1])  # x^3 - 3x^2 + 2x

    def test_scale(self):
        p = Polynomial([FieldElement(1, self.field), FieldElement(2, self.field), FieldElement(1, self.field)], self.field)  # 1 + 2x + x^2
        scaled = p.scale(2)
        self.assertEqual(scaled.poly.tolist(), [1, 4, 4])  # 1 + 4x + 4x^2

    def test_colinearity(self):
        points = [(FieldElement(0, self.field), FieldElement(1, self.field)),
                  (FieldElement(1, self.field), FieldElement(2, self.field)),
                  (FieldElement(2, self.field), FieldElement(3, self.field))]
        self.assertTrue(Polynomial.test_colinearity(points, self.field))

        non_colinear_points = [(FieldElement(0, self.field), FieldElement(1, self.field)),
                               (FieldElement(1, self.field), FieldElement(2, self.field)),
                               (FieldElement(2, self.field), FieldElement(5, self.field))]
        self.assertFalse(Polynomial.test_colinearity(non_colinear_points, self.field))

if __name__ == '__main__':
    unittest.main()