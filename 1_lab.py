import matplotlib.pyplot as plt
import numpy as np


def step(right, left, number_of_dots):
    return (right - left) / number_of_dots


def main_equation(x):
    return np.sqrt(2*x) - np.cos(x/3)


def f(x):
    return (np.sin(x/3)/3) + (np.sqrt(2) / (2 * np.sqrt(x)))


def ff(x):
    return (np.cos(x/3)/9) + (np.sqrt(2) / (4 * np.sqrt(x**3)))


def main_calculation():
    stepped_dot = 0.5
    dots_with_step = []
    for i in range(n + 1):
        dots_with_step.append(stepped_dot)
        stepped_dot += step(right_edge, left_edge, n)
    dots_of_main_equation = [main_equation(i) for i in range(10)]
    dots_of_stepped_equation = [main_equation(i) for i in dots_with_step]
    plt.plot(dots_of_stepped_equation)
    print('Dots of main equation')
    print(dots_of_main_equation)
    print('Dots of stepped equation')
    print(dots_of_stepped_equation)
    plt.show()


def tangent_method():
    def newton(xk):
        return xk - main_equation(xk) / f(xk)
    x0 = 0.5
    values = []
    while np.abs(newton(x0) - x0) >= 10**(-7):
        xk1 = x0
        x0 = newton(xk1)
        values.append(round(x0, 7))

    print('Tangment method')
    print(values)
    print('Error:', values[-1] - values[-2])


def chord_method():
    x0 = left_edge
    values = []

    def chord_formula(xk):
        return xk - ((main_equation(xk)) / (main_equation(right_edge) - main_equation(xk)) * (right_edge - xk))

    while np.abs(chord_formula(x0) - x0) >= 10**(-7):
        xk1 = x0
        x0 = chord_formula(xk1)
        values.append(round(x0, 7))

    print('Chord method')
    print(values)
    print('Error:', values[-2] - values[-1], '/n/n')


def secant_method():
    x0 = 0.1
    x1 = 0.4
    values = [x0, x1]

    def secant_formula(xk, xk_1):
        c = (xk - xk_1)
        b = (main_equation(xk) - main_equation(xk_1))
        a = xk - (main_equation(xk) / 1.0 * b) * c
        return a

    i = 1
    while np.abs(secant_formula(values[i], values[i-1]) - values[i-1]) >= 10**(-7):
        xk1 = values[i]
        x0 = secant_formula(xk1, values[i-1])
        values.append(round(x0, 7))
        i += 1

    print('Secant method')
    print(values)
    print('Error:', values[-2] - values[-3])


def finite_difference_method():
    def finite_difference_formula(xk, h):
        return xk - (h * main_equation(xk)) / (main_equation(xk+h) - main_equation(xk))

    x0 = 0.1
    h0 = 0.05
    values = []
    while np.abs(finite_difference_formula(x0, h0) - x0) >= 10**(-7):
        xk1 = x0
        x0 = finite_difference_formula(xk1, h0)
        values.append(round(x0, 7))

    print('Finite difference method')
    print(values)
    print('Error:', values[-1] - values[-2])


def steffans_method():
    def steffans_formula(xk):
        return xk - ((main_equation(xk))**2) / (main_equation(xk+main_equation(xk)) - main_equation(xk))

    values = []
    x0 = 0.4
    while np.abs(steffans_formula(x0) - x0) >= 10**(-7):
        xk1 = x0
        x0 = steffans_formula(xk1)
        values.append(round(x0, 7))

    print('Steffans method')
    print(values)
    print('Error:', values[-1] - values[-2])


def simple_iterations_method():
    def simple_iterations_formula(xk, t):
        return xk - t * main_equation(xk)

    t0 = 0.5
    x0 = 0.4
    values = []
    while np.abs(simple_iterations_formula(x0, t0) - x0) >= 10**(-7):
        xk1 = x0
        x0 = simple_iterations_formula(xk1, t0)
        values.append(round(x0, 7))

    print('Simple iterations method')
    print(values)
    print('Error:', values[-1] - values[-2])


left_edge = 0.1
right_edge = 0.5
n = 10

tangent_method()
chord_method()
secant_method()
finite_difference_method()
steffans_method()
simple_iterations_method()
main_calculation()
