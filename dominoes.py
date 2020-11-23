from math import sin, cos, sqrt, asin, acos, atan, pi
import matplotlib.pyplot as plt
import numpy as np


# mgh / I = 3gh / (h^2 + d^2) = 3g/h * 1/(1 + gamma^2)


def previous_angle(angle, beta, gamma):
    """
    theta_{i-1}(theta_i)
    :param angle: theta_i
    :return: theta_{i-1}
    """
    return angle - asin(beta * sin(angle) - gamma)


def xi(theta, beta, gamma):
    return 1 - beta * cos(theta) / sqrt(1 - (beta * sin(theta) - gamma) ** 2)


def sum_of_sin(N, theta_N, beta, gamma):
    """
    sum from i=1 to N sin(theta_i + atan(gamma))
    :param N: number of involved dominoes
    :param theta_N: angle of the last domino
    :param beta:
    :param gamma:
    :return:
    """
    summ = 0
    theta_i = theta_N
    for i in range(N, 0, -1):
        summ += sin(theta_i + atan(gamma))
        theta_i = previous_angle(angle=theta_i, beta=beta, gamma=gamma)
    return summ


def sum_of_squared_xi(N, theta_N, beta, gamma):
    """
    sum from i=1 to N xi(theta_N)^2 * xi(theta_{N-1})**2 * ... * xi(theta_{i+1})**2
    :param N: number of involved dominoes
    :param theta_N: angle of the last domino
    :param beta:
    :param gamma:
    :return:
    """
    summ = 0
    prod = 1
    theta_i = theta_N
    for i in range(N, 0, -1):
        summ += prod
        prod *= xi(theta=theta_i, beta=beta, gamma=gamma) ** 2
        theta_i = previous_angle(angle=theta_i, beta=beta, gamma=gamma)
    return summ


def omega_intermediate(theta, omega_0, N, beta, gamma, h):
    """
    :param theta: theta_N
    :param omega_0: \dot{theta_N}(theta_N = pi/2)
    :param N: domino number
    :return: \dot{theta_N}
    """
    try:
        return sqrt(
            (3 * g / h * 1 / sqrt(1 + gamma ** 2) *
             (sum_of_sin(N=N, theta_N=pi / 2, beta=beta, gamma=gamma) - sum_of_sin(N=N, theta_N=theta, beta=beta,
                                                                                   gamma=gamma))
             + omega_0 ** 2 * sum_of_squared_xi(N=N, theta_N=pi / 2, beta=beta, gamma=gamma))
            / sum_of_squared_xi(N=N, theta_N=theta, beta=beta, gamma=gamma))
    except ValueError as e:
        print("WARNING: small energy; answers are incorrect")
        print('It is less than zero ', (3 * g / h * 1 / sqrt(1 + gamma ** 2) *
             (sum_of_sin(N=N, theta_N=pi / 2, beta=beta, gamma=gamma) - sum_of_sin(N=N, theta_N=theta, beta=beta,
                                                                                   gamma=gamma))
             + omega_0 ** 2 * sum_of_squared_xi(N=N, theta_N=pi / 2, beta=beta, gamma=gamma))
            / sum_of_squared_xi(N=N, theta_N=theta, beta=beta, gamma=gamma))
        print('sum_of_squared_xi ', sum_of_squared_xi(N=N, theta_N=pi / 2, beta=beta, gamma=gamma), 'N', N)
        return 1


def omega_after_collision(omega_before_collision, N, beta, gamma):
    """
    calculate collision
    :param omega_before_collision:
    :param N: number domino that pushes next one
    :param beta:
    :param gamma:
    :return:
    """
    return omega_before_collision * sum_of_squared_xi(N=N, theta_N=previous_angle(angle=pi/2, beta=beta, gamma=gamma),
                                                      beta=beta, gamma=gamma) / sum_of_squared_xi(N=N+1,
                                                                                                  theta_N=pi / 2,
                                                                                                  beta=beta,
                                                                                                  gamma=gamma)


def asymptotic_velocity(n, omega_initial, beta, gamma):
    times_min = []
    times_max = []
    times = [0]
    # t_0 = float(input("t_0 = "))
    omega_last_in_vertical_position = omega_initial
    for i in range(1, n + 1):  # i - number of involved dominoes
        omegas = [omega_last_in_vertical_position]
        t_min, t_max = 0, 0
        for j in range(1, DEPTH + 1):
            omega_last = omega_intermediate(theta=pi / 2 - j / DEPTH * delta, omega_0=omega_last_in_vertical_position,
                                            N=i, beta=beta, gamma=gamma, h=h)
            omegas.append(omega_last)
        for j in range(len(omegas) - 1):
            t_min += d_theta / omegas[j + 1]
            t_max += d_theta / omegas[j]
        times_min.append(t_min)
        times_max.append(t_max)
        times.append(times[-1] + (t_min + t_max) / 2)
        omega_last_in_vertical_position = omega_after_collision(omega_before_collision=omegas[-1], N=i, beta=beta,
                                                                gamma=gamma)

    # print(*map(lambda x: x + t_0 - times[1], times[2:]), sep='\n')
    v_min, v_max = s / times_min[-1], s / times_max[-1]  # v_max, v_min
    # return v_max, v_min
    return v_max, v_min



g = 9.8
n = 20
DEPTH = 100

h = 0.07
# s = 0.04
d = 0.01
# beta = s / h
gamma = d / h
omega_initial = sqrt(3 * g / h * 1 / (1 + gamma ** 2) * (sqrt(1 + gamma ** 2) - 1)) + 0.01
# omega_initial += 5.5
# print(omega_initial)
# print()


# delta = pi/2 - previous_angle(angle=pi/2, beta=beta, gamma=gamma)
# d_theta = delta / DEPTH


def main():
    global d_theta, delta, s
    s_array = np.arange(d + h * sin(2 * atan(gamma)) + 0.001, h + 0.001, 0.001)
    v_max_array = []
    v_min_array = []
    counter = 0
    for s in s_array:
        counter += 1
        print((counter * 100) // len(s_array), '%', sep=' ')

        beta = s / h
        delta = pi / 2 - previous_angle(angle=pi / 2, beta=beta, gamma=gamma)
        d_theta = delta / DEPTH

        v_max, v_min = asymptotic_velocity(n=n, omega_initial=omega_initial, beta=beta, gamma=gamma)
        v_max_array.append(v_max)
        v_min_array.append(v_min)
    plt.plot(s_array, v_max_array, s_array, v_min_array)
    plt.grid()
    plt.show()


def main1():
    global d_theta, delta, s
    s = 0.06
    beta = s / h
    delta = pi / 2 - previous_angle(angle=pi / 2, beta=beta, gamma=gamma)
    d_theta = delta / DEPTH

    print(asymptotic_velocity(n=n, omega_initial=omega_initial, beta=beta, gamma=gamma))


if __name__ == '__main__':
    main1()

# for i in range(1, n + 1):  # i - number of involved dominoes
#     omegas = [omega_last_in_vertical_position]
#     t_min, t_max = 0, 0
#     for j in range(1, DEPTH + 1):
#         omega_last = omega_intermediate(theta=pi/2 - j / DEPTH * delta, omega_0=omega_last_in_vertical_position, N=i, beta=beta, gamma=gamma, h=h)
#         omegas.append(omega_last)
#     for j in range(len(omegas) - 1):
#         t_min += d_theta / omegas[j + 1]
#         t_max += d_theta / omegas[j]
#     times_min.append(t_min)
#     times_max.append(t_max)
#     omega_last_in_vertical_position = omega_after_collision(omega_before_collision=omegas[-1], N=i, beta=beta, gamma=gamma)
# print(times_min)
# print(times_max)
# print("v_max = ", s / times_min[-1])
# print("v_min = ", s / times_max[-1])
