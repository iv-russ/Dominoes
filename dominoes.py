from math import sin, cos, tan, sqrt, asin, acos, atan, pi
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


def zeta(theta, beta, gamma):
    constant = beta * sin(theta) - gamma
    return beta * (sin(theta) * (1 - constant ** 2) - beta * cos(theta) ** 2 * constant) / (1 - constant ** 2) ** (
            3 / 2)


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
    sum from i=1 to N xi(theta_N)^2 * xi(theta_{N-1})^2 * ... * xi(theta_{i+1})^2
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


def A_eff(theta_k, beta, gamma, mu):
    theta_k_minus_1 = previous_angle(angle=theta_k, beta=beta, gamma=gamma)
    return sqrt(1 + (beta - gamma / sin(theta_k)) ** 2 - 2 * (beta - gamma / sin(theta_k)) * cos(theta_k_minus_1)) - \
           gamma / tan(theta_k) - mu * gamma


def B_eff(theta_k, beta, gamma, mu):
    return sqrt(1 - (beta * sin(theta_k) - gamma) ** 2) + mu * (beta * sin(theta_k) - gamma)


def G(theta_k, gamma):
    return sqrt(1 + gamma ** 2) / 2 * cos(theta_k + atan(gamma))


def sum_of_xi_and_coefs(N, theta_N, beta, gamma, mu):
    summ = 0
    prod = 1
    theta_i = theta_N
    for i in range(N, 0, -1):
        summ += prod
        prod *= xi(theta=theta_i, beta=beta, gamma=gamma) * \
                A_eff(theta_k=theta_i, beta=beta, gamma=gamma, mu=mu) / \
                B_eff(theta_k=theta_i, beta=beta, gamma=gamma, mu=mu)
        theta_i = previous_angle(angle=theta_i, beta=beta, gamma=gamma)
    return summ


def sum_of_G(N, theta_N, beta, gamma, mu):
    summ = 0
    prod = 1
    theta_i = theta_N
    for i in range(N, 0, -1):
        summ += G(theta_k=theta_i, gamma=gamma) * prod
        prod *= A_eff(theta_k=theta_i, beta=beta, gamma=gamma, mu=mu) / \
                B_eff(theta_k=theta_i, beta=beta, gamma=gamma, mu=mu)
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
    except ValueError:
        print("WARNING: small energy; answers are incorrect")
        print('It is less than zero ', (3 * g / h * 1 / sqrt(1 + gamma ** 2) *
                                        (sum_of_sin(N=N, theta_N=pi / 2, beta=beta, gamma=gamma) - sum_of_sin(N=N,
                                                                                                              theta_N=theta,
                                                                                                              beta=beta,
                                                                                                              gamma=gamma))
                                        + omega_0 ** 2 * sum_of_squared_xi(N=N, theta_N=pi / 2, beta=beta, gamma=gamma))
              / sum_of_squared_xi(N=N, theta_N=theta, beta=beta, gamma=gamma))
        return 1


def omega_after_collision(omega_before_collision, N, beta, gamma, mu):
    """
    calculate collision
    :param omega_before_collision:
    :param N: number domino that pushes next one
    :param beta:
    :param gamma:
    :return:
    """
    # print(sum_of_xi_and_coefs(N=N, theta_N=previous_angle(angle=pi / 2, beta=beta, gamma=gamma), beta=beta, gamma=gamma, mu=mu))
    # print(sum_of_xi_and_coefs(N=N + 1, theta_N=pi / 2, beta=beta, gamma=gamma, mu=mu))
    # print(sum_of_xi_and_coefs(N=N, theta_N=previous_angle(angle=pi / 2, beta=beta, gamma=gamma), beta=beta, gamma=gamma,
    #                           mu=0))
    # print(sum_of_xi_and_coefs(N=N + 1, theta_N=pi / 2, beta=beta, gamma=gamma, mu=0))
    print(omega_before_collision)
    omega_after = omega_before_collision * A_eff(theta_k=pi / 2, beta=beta, gamma=gamma, mu=mu) / B_eff(theta_k=pi / 2,
                                                                                                 beta=beta, gamma=gamma,
                                                                                                 mu=mu) * \
           sum_of_xi_and_coefs(N=N, theta_N=previous_angle(angle=pi / 2, beta=beta, gamma=gamma), beta=beta,
                               gamma=gamma, mu=mu) / sum_of_xi_and_coefs(N=N + 1,
                                                                         theta_N=pi / 2,
                                                                         beta=beta,
                                                                         gamma=gamma,
                                                                         mu=mu)
    print(omega_after)
    return omega_before_collision * A_eff(theta_k=pi / 2, beta=beta, gamma=gamma, mu=mu) / B_eff(theta_k=pi / 2,
                                                                                                 beta=beta, gamma=gamma,
                                                                                                 mu=mu) * \
           sum_of_xi_and_coefs(N=N, theta_N=previous_angle(angle=pi / 2, beta=beta, gamma=gamma), beta=beta,
                               gamma=gamma, mu=mu) / sum_of_xi_and_coefs(N=N + 1,
                                                                         theta_N=pi / 2,
                                                                         beta=beta,
                                                                         gamma=gamma,
                                                                         mu=mu)


def time_of_falling_and_omega(omega_0, N, beta, gamma, h, dt):
    t = 0
    theta_1 = previous_angle(angle=pi / 2, beta=beta, gamma=gamma)
    thetas = [pi / 2]
    for _ in range(N - 1):
        thetas.append(previous_angle(angle=thetas[-1], beta=beta, gamma=gamma))
    omegas = [omega_0]
    for i in range(1, N):
        omegas.append(xi(theta=thetas[i - 1], beta=beta, gamma=gamma) * omegas[-1])
    while thetas[0] >= theta_1:
        # sum_of_G = 0
        # prod = 1
        # for i in range(N):
        #     sum_of_G += G(theta_k=thetas[i], gamma=gamma) * prod
        #     prod *= A_eff(theta_k=thetas[i], beta=beta, gamma=gamma, mu=mu)

        s_G = sum_of_G(N=N, theta_N=thetas[0], beta=beta, gamma=gamma, mu=mu)

        ddot_theta = g / h * 3 / (1 + gamma ** 2) * s_G / sum_of_xi_and_coefs(N=N, beta=beta, gamma=gamma,
                                                                              mu=mu,
                                                                              theta_N=thetas[0])
        thetas[0] -= omegas[0] * dt
        for i in range(1, N):
            thetas[i] = previous_angle(angle=thetas[i - 1], beta=beta, gamma=gamma)
        omegas[0] += ddot_theta * dt
        for i in range(1, N):
            omegas[i] = xi(theta=thetas[i - 1], beta=beta, gamma=gamma) * omegas[i - 1]

        t += dt
    return t, omegas[0]


def asymptotic_velocity(n, omega_initial, beta, gamma, mu):
    t = 0
    omega_last_in_vertical_position = omega_initial
    for i in range(1, n + 1):  # i - number of involved dominoes
        # omegas = [omega_last_in_vertical_position]
        # t_min, t_max = 0, 0
        # for j in range(1, DEPTH + 1):
        #     omega_last = omega_intermediate(theta=pi / 2 - j / DEPTH * delta, omega_0=omega_last_in_vertical_position,
        #                                     N=i, beta=beta, gamma=gamma, h=h)
        #     omegas.append(omega_last)
        # for j in range(len(omegas) - 1):
        #     t_min += d_theta / omegas[j + 1]
        #     t_max += d_theta / omegas[j]
        # times_min.append(t_min)
        # times_max.append(t_max)
        # times.append(times[-1] + (t_min + t_max) / 2)
        dt = 0.0001
        t, omega = time_of_falling_and_omega(omega_0=omega_last_in_vertical_position, N=i, beta=beta, gamma=gamma, h=h,
                                             dt=dt)
        omega_last_in_vertical_position = omega_after_collision(omega_before_collision=omega, N=i, beta=beta,
                                                                gamma=gamma, mu=mu)

    # print(*map(lambda x: x + t_0 - times[1], times[2:]), sep='\n')
    s = beta * h
    v = s / t  # velocity
    return v


g = 9.8
mu = 0
n = 40
DEPTH = 100

h = 0.117
# s = 0.04
d = 0.008
# beta = s / h
gamma = d / h
omega_initial = sqrt(3 * g / h * 1 / (1 + gamma ** 2) * (sqrt(1 + gamma ** 2) - 1)) + 0.01


def main():
    global d_theta, delta
    s_array = np.arange(d + h * sin(1.5 * atan(gamma)) + 0.001, h + 0.001, 0.001)
    # s_array = np.arange(0.02, 0.10, 0.001)
    v_array = []
    counter = 0
    for s in s_array:
        counter += 1
        print((counter * 100) // len(s_array), '%', sep=' ')

        beta = s / h
        delta = pi / 2 - previous_angle(angle=pi / 2, beta=beta, gamma=gamma)
        d_theta = delta / DEPTH

        v = asymptotic_velocity(n=n, omega_initial=omega_initial, beta=beta, gamma=gamma, mu=mu)
        v_array.append(v)
    plt.plot(s_array, v_array)
    plt.grid()
    plt.show()


def main1():
    global d_theta, delta, s
    s = 0.02
    beta = s / h
    delta = pi / 2 - previous_angle(angle=pi / 2, beta=beta, gamma=gamma)
    d_theta = delta / DEPTH
    print(asymptotic_velocity(n=n, omega_initial=omega_initial, beta=beta, gamma=gamma, mu=mu))

    # s_array = np.arange(0.02, 0.10, 0.001)
    # v_array = []
    # for s in s_array:
    #     beta = s/h
    #     v_array.append(asymptotic_velocity(n=n, omega_initial=omega_initial, beta=beta, gamma=gamma, mu=mu))
    # plt.plot(s_array, v_array)
    # plt.grid()
    # plt.show()
    # print(asymptotic_velocity(n=n, omega_initial=omega_initial, beta=beta, gamma=gamma, mu=mu))


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
