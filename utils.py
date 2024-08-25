from tinygrad import Tensor

M31 = 2**31 - 1
M31SQ = (2**31 - 1) ** 2
HALF = 2**30


def array(x):
    return Tensor(x)


def zeros(x):
    return Tensor.zeros(x)


def tobytes(x):
    return x.numpy().tobytes()


def arange(*args):
    return Tensor.arange(*args)


def append(*args):
    return args[0].cat(args[1:])


def log2(x):
    assert x & (x - 1) == 0
    return x.bit_length() - 1


def pad_to(arr, new_len):
    padding = zeros((new_len - arr.shape[0],) + arr.shape[1:])
    return append(arr, padding)


def confirm_max_degree(coeffs, bound):
    return not coeffs[bound:].any()
