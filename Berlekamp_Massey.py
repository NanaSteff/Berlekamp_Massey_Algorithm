import sympy as sy              # importation of sympy
sy.init_printing()              # output latex display
x = sy.Symbol('x')


def Berlekamp(seq):
    f = 1
    deg_f_list = []
    n = len(seq)
    f_list = [1]
    for i in range(n):
        d = 0
        print("step i = ", i)
        deg_f = sy.degree(f)
        deg_f_list.append(deg_f)
        coeffs_f = sy.Poly(f, x).all_coeffs()   
        if coeffs_f != []:
            coeffs_f = coeffs_f[-1::-1]
        for j in range(deg_f+1):
            d = d+coeffs_f[deg_f-j] * seq[i-j]
            d = d % 2
        if d == 0:
            f = sy.Poly(f, x).trunc(2)
            f_list.append(f)
        else:
            L_index = [k for k in range(
                len(deg_f_list)-1) if deg_f_list[k] < deg_f_list[k+1]]
            if len(L_index) == 0:
                m = -1
                Lm = 0
                fm = 1
            else:
                m = max(L_index)
                Lm = deg_f_list[m]
                fm = f_list[m]
            Li = deg_f_list[i]
            if (m-Lm) >= (i-Li):
                f = f+x**((m-Lm)-(i-Li))*fm
                f = sy.Poly(f, x).trunc(2)
                f_list.append(f)
            else:
                f = (x**((i-Li)-(m-Lm)))*(f)+fm
                f = sy.Poly(f, x).trunc(2)
                f_list.append(f)
        print('f({})='.format(i+1), f)

    return f


# driver code
if __name__ == "__main__":
    seq = [1, 1, 0, 1, 0, 0, 1]
    Berlekamp(seq)