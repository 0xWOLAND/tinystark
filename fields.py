from tinygrad import Tensor
import numpy as np

class FieldElement: 
    def __init__(self, value, field):
        self.value = value
        self.field: Field = field

    def __add__(self, right):
        return self.field.add(self, right)
    
    def __sub__(self, right):
        return self.field.sub(self, right)
    
    def __mul__(self, right):
        return self.field.mul(self, right)
    
    def __truediv__(self, right):
        return self.field.divide(self, right)
    
    def __neg__(self):
        return self.field.negate(self)
    
    def inverse(self):
        return self.field.inverse(self)
    
    def __eq__(self, right):
        return self.value == right.value
    
    def __ne__(self, right):
        return self.value != right.value
    
    def __str__(self):
        return str(self.value)
    
    def __bytes__(self):
        return bytes(str(self).encode())
    
    def __xor__(self, exp):
        acc = FieldElement(1, self.field)
        val = FieldElement(self.value, self.field)
        for i in reversed(range(len(bin(exp)[2:]))):
            acc = acc * acc
            if (1 << i) & exp != 0:
                acc * val
        return acc

    def is_zero(self):
        return self.value == 0

def xgcd(x, y):
    r, rr = x, y
    s, ss = 1, 0
    t, tt = 0, 1

    while rr != 0:
        q = r // rr
        r, rr = rr, r - q * rr
        s, ss = ss, s - q * ss
        t, tt = tt, t - q * tt

    return s, t, r


class Field:
    def __init__(self, p):
        self.p = p

    def zero(self):
        return FieldElement(0, self)
    
    def one(self):
        return FieldElement(1, self)
    
    def add(self, left, right):
        return FieldElement((left.value + right.value) % self.p, self)
    
    def sub(self, left, right):
        return FieldElement((left.value - right.value) % self.p, self)
    
    def mul(self, left, right):
        return FieldElement((left.value * right.value) % self.p, self)
    
    def divide(self, left, right):
        if right.is_zero():
            raise ZeroDivisionError('Cannot divide by zero')
        a, b, g = xgcd(right.value, self.p)
        return FieldElement((left.value * a) % self.p, self)

    def inverse(self, operand):
        a, b, g = xgcd(operand.value, self.p)
        return FieldElement(a, self)

    def negate(self, operand):
        return FieldElement((self.p -operand.value) % self.p, self)
    
    def generator( self ):
        assert(self.p == 1 + 407 * ( 1 << 119 )), "Do not know generator for other fields beyond 1+407*2^119"
        return FieldElement(85408008396924667383611388730472331217, self)
        
    def primitive_nth_root( self, n ):
        if self.p == 1 + 407 * ( 1 << 119 ):
            assert(n <= 1 << 119 and (n & (n-1)) == 0), "Field does not have nth root of unity where n > 2^119 or not power of two."
            root = FieldElement(85408008396924667383611388730472331217, self)
            order = 1 << 119
            while order != n:
                root = root^2
                order = order/2
            return root
        else:
            assert(False), "Unknown field, can't return root of unity."

    def sample( self, byte_array ):
        acc = 0
        for b in byte_array:
            acc = (acc << 8) ^ int(b)
        return FieldElement(acc % self.p, self)
    
    def to_field_element( self, value ):
        return FieldElement(value % self.p, self)

class Polynomial:
    def __init__(self, poly, field: 'Field'):
        if isinstance(poly, list):
            if all(isinstance(coeff, (int, float)) for coeff in poly):
                self.poly = Tensor(poly)
            elif all(isinstance(coeff, FieldElement) for coeff in poly):
                self.poly = Tensor([x.value for x in poly])
        else:
            self.poly = poly
        self.field = field
    
    def degree(self):
        if self.poly.shape[0] == 0 or np.all(self.poly.numpy() == 0):
            return -1

        for i, coeff in enumerate(self.poly.tolist()):
            if coeff != 0:
                return i
        return 0

    def __add__(self, right: 'Polynomial') -> 'Polynomial':
        if self.degree() == -1:
            return right
        if right.degree() == -1:
            return self
    
        n = max(self.degree(), right.degree())
        lhs = self.poly.pad(((0, n),))
        rhs = right.poly.pad(((0, n),))

        tot = (lhs + rhs).numpy().astype(int)
        tot = np.vectorize(self.field.to_field_element)(tot)
        return Polynomial(tot, self.field)

    def __sub__(self, right: 'Polynomial') -> 'Polynomial':
        if self.degree() == -1:
            return right
        if right.degree() == -1:
            return self

        n = max(self.degree(), right.degree())
        lhs = self.poly.pad(((0, n),))
        rhs = right.poly.pad(((0, n),))

        diff = (lhs - rhs).numpy().astype(int)
        diff = np.vectorize(self.field.to_field_element)(diff)
        return Polynomial(diff, self.field)

    def __mul__(self, other):
        if self.poly.shape[0] == 0 or other.poly.shape[0] == 0:
            return Polynomial([], self.field)
        
        # Create a 2D tensor for efficient multiplication
        m, n = self.poly.shape[0], other.poly.shape[0]
        result_size = m + n - 1
        
        # Create a matrix where each row is a shifted version of other.poly
        matrix = Tensor.zeros((m, result_size), dtype=self.poly.dtype)
        for i in range(m):
            matrix[i, i:i+n] = other.poly.reshape(n)  # Ensure other.poly is 1D
        
        # Perform matrix multiplication
        result = Tensor.matmul(self.poly.reshape(1, -1), matrix).reshape(-1)
        
        return Polynomial(result, self.field)
    
    def __truediv__(self, right: 'Polynomial') -> 'Polynomial':
        quo, rem = Polynomial.divide(self, right, self.field)
        assert(rem.is_zero()), "Division is not exact"
        return quo

    def __mod__(self, right: 'Polynomial') -> 'Polynomial':
        quo, rem = Polynomial.divide(self, right)
        return rem

    def __xor__(self, exp):
        if self.is_zero():
            return Polynomial([])
        if exp == 0:
            return Polynomial([self.field.one()])
        acc = Polynomial([self.field.one()])
        for i in reversed(range(len(bin(exp)[2:]))):
            acc = acc * acc
            if (1 << i) & exp != 0:
                acc = acc * self
        return acc

    def __neg__(self):
        self.poly = self.poly.neg() 

    def __eq__(self, right: 'Polynomial') -> bool:
        return self.poly == right.poly

    def __ne__(self, right: 'Polynomial') -> bool:
        return self.poly != right.poly
    
    def is_zero(self):
        return self.poly == Tensor.zeros(self.poly.shape)

    def leading_coefficient(self):
        return self.poly[self.degree()]

    def divide(numerator, denominator, field):
        if denominator.degree() == -1:
            return None
        if numerator.degree() < denominator.degree():
            return (Polynomial([]), numerator)

        n_degree, d_degree = numerator.degree(), denominator.degree()
        quotient_degree = n_degree - d_degree + 1

        # Convert polynomials to tensors for efficient operations
        n_coeffs = numerator.poly
        d_coeffs = denominator.poly
        
        # Prepare tensors for quotient and remainder
        quotient_coeffs = Tensor.zeros(quotient_degree, dtype=n_coeffs.dtype)
        remainder_coeffs = n_coeffs

        # Precompute inverse of leading coefficient of denominator
        d_leading_inv = 1 / d_coeffs[-1]

        for i in range(quotient_degree):
            if remainder_coeffs.shape[0] < d_coeffs.shape[0]:
                break
            
            # Compute the current quotient coefficient
            coeff = remainder_coeffs[-1] * d_leading_inv
            quotient_coeffs[-(i+1)] = coeff

            # Compute the polynomial to subtract
            shift = remainder_coeffs.shape[0] - d_coeffs.shape[0]
            subtractee = Tensor.zeros(shift, dtype=n_coeffs.dtype).cat(d_coeffs * coeff)
            
            # Perform the subtraction
            remainder_coeffs = remainder_coeffs - subtractee

            # Trim trailing zeros from remainder
            remainder_coeffs = remainder_coeffs[:(remainder_coeffs != 0).sum().item()]

        quotient = Polynomial(quotient_coeffs, field)
        remainder = Polynomial(remainder_coeffs, field)

        return quotient, remainder
    
    def evaluate(self, point: FieldElement) -> FieldElement:
        if len(self.poly) == 0:
            return self.field.zero()

        # Create a tensor of powers of the point
        powers = Tensor.arange(len(self.poly), dtype=self.poly.dtype)
        point_powers = point.value ** powers

        # Perform element-wise multiplication and sum
        result = Tensor.sum(self.poly * point_powers).item()

        return FieldElement(result, self.field)

    def evaluate_domain(self, domain):
        if len(self.poly) == 0:
            return [self.field.zero() for _ in domain]

        # Convert poly and domain to tensors
        coeffs = Tensor(self.poly)
        domain_tensor = Tensor(domain)

        # Create a 2D tensor of powers of domain elements
        degrees = Tensor.arange(len(coeffs), dtype=coeffs.dtype)
        domain_powers = domain_tensor.unsqueeze(1) ** degrees.unsqueeze(0)

        # Perform batch multiplication and sum
        results = Tensor.matmul(domain_powers, coeffs.unsqueeze(1)).squeeze(1)

        return results.tolist()

    def interpolate_domain(domain, values, field):
        assert len(domain) == len(values), "number of elements in domain does not match number of values -- cannot interpolate"
        assert len(domain) > 0, "cannot interpolate between zero points"
        
        n = len(domain)
        
        # Convert domain and values to tensors
        domain_tensor = Tensor([d.value for d in domain])
        values_tensor = Tensor([v.value for v in values])
        
        # Create a tensor for (x - x_j) terms
        x_minus_xj = domain_tensor.reshape(1, n) - domain_tensor.reshape(n, 1)
        
        # Replace diagonal elements with ones to avoid division by zero
        x_minus_xj = x_minus_xj + Tensor.eye(n)
        
        # Compute the inverse of (x_i - x_j) for i != j
        inv_diff = x_minus_xj.reciprocal()
        inv_diff = inv_diff * (1 - Tensor.eye(n))
        
        # Compute the product of all (x_i - x_j) for i != j
        prod_diff = inv_diff.prod(axis=1)
        
        # Compute the poly of the interpolation polynomial
        coeff = values_tensor * prod_diff
        
        # Generate powers of x for each term
        x_powers = Tensor.arange(n).reshape(1, n).repeat((n, 1))
        
        # Compute the final poly
        final_coeff = Tensor.zeros(n)
        for i in range(n):
            term = coeff[i] * ((-domain_tensor[i]) ** x_powers[i])
            
            # Perform polynomial multiplication without conv1d
            expanded_term = Tensor.zeros((n, 2*n-1))
            for j in range(n):
                expanded_term[i, j:j+n] += term[j] * inv_diff[i]
            
            final_coeff += expanded_term.sum(axis=0)[:n]
        
        # Convert final poly back to field elements
        final_coeff_field = [field.to_field_element(int(c)) for c in final_coeff.numpy()]
        
        return Polynomial(final_coeff_field)
        
    def zeroifier_domain(domain, field):
        x = Polynomial([field.zero(), field.one()], field)
        acc = Polynomial([field.one()], field)
        for d in domain:
            acc = acc * (x - Polynomial([d], field))
        return acc
    
    import numpy as np

    def scale(self, scaling_factor):
        if len(self.poly) == 0:
            return Polynomial([], self.field)
        
        # Convert scaling_factor to a field element if it's not already
        if not isinstance(scaling_factor, FieldElement):
            scaling_factor = self.field.to_field_element(scaling_factor)
        
        # Create a tensor of powers of the scaling_factor
        exponents = Tensor.arange(len(self.poly))
        scaling_factors = Tensor([scaling_factor.value]).pow(exponents)
        
        # Multiply the poly by the scaling factors
        scaled_poly_tensor = self.poly * scaling_factors
        
        # Convert the scaled poly back to field elements
        scaled_poly = [self.field.to_field_element(int(c)) for c in scaled_poly_tensor.numpy()]
        
        return Polynomial(scaled_poly, self.field)

    def test_colinearity(points, field):
        domain = [p[0] for p in points]
        values = [p[1] for p in points]
        polynomial = Polynomial.interpolate_domain(domain, values, field)
        return polynomial.degree() <= 1



def main():
    p = 1 + 407 * ( 1 << 119 ) # 1 + 11 * 37 * 2^119
    return Field(p)

if __name__ == '__main__':
    main()