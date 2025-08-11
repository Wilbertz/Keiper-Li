from mpmath import *

def get_coefficients(filename):
    with open(filename) as csvfile:
        lines = csvfile.readlines()
        result = []
        for line in lines:
            number = mpf(line.split(';')[1])
            result.append(number)
        return result

def get_keiper_li():
    return get_coefficients("Keiper-Li-1-30.csv")

def get_stieltjes():
    return get_coefficients("Stieltjes-1-30.csv")

if __name__ == "__main__":
    mp.prec = 100
    mp.dps = 100
    coefficients = get_stieltjes()
    print(coefficients[3])
