
    """

    R_ = math.sqrt(V_.x_**2 + V_.y_**2 + V_.z_**2)

    # Функция - подинтегральное выражение
    def f_(t_1):
        tau_ = t_1
        return math.exp((-v_**2 * tau_)/(4*a_) - (R_**2)/(4*a_*tau_))*(1/(tau_**(3/2)))

    i_ = scipy.integrate.quad(f_, 0, t_)
    return T_n + (2*q_)/(cp_*math.sqrt((4*math.pi*a_)**3)) * math.exp((-v_*V_.x_)/(2*a_)) * i_[0]

# Формула для состояния температурного поля при воздействии быстро движущегося линейного источника
def T_2(T_n, V_, t_, q_, cp_, a_, v_, lambda_, delta_):
    """
    T_n - начальная температура изделия
    x_ - координата x
    t_ - время
    q_ - мощность
    cp_ - теплоемкость материала
    a_ - коэффициент температуропроводности
    v_ - скорость сварки
    lambda_ - коэффициент теплопроводности
    delta_ - толщина
    """

    # Функция - подинтегральное выражение
    def f_(t_1):
        tau_ = t_1
        return math.exp( (-v_**2 * tau_)/(4*a_) - (2*lambda_*tau_)/(cp_*delta_) - (V_.x_**2+V_.y_**2)/(4*a_*tau_) ) * (1/tau_)

    i_ = scipy.integrate.quad(f_, 0, t_, limit=10)
    return T_n + (q_)/(4*math.pi*lambda_*delta_) * math.exp((-v_*V_.x_)/(2*a_)) * i_[0]

# Функция ψ(x, y, z, v, t, q)
def PSI_xyzvtq(T_n, V_, t_, q_, cp_, a_, v_, lambda_, delta_):
    """
    T_n - начальная температура изделия
    x_ - координата x
    t_ - время
    q_ - мощность
    cp_ - теплоемкость материала
    a_ - коэффициент температуропроводности
    v_ - скорость сварки
    lambda_ - коэффициент теплопроводности
    delta_ - толщина
    T_t - температура в стадии теплонасыщения
    """

    T_1_result = T_1(T_n, V_, t_, q_, cp_, a_, v_)
    T_2_result = T_2(T_n, V_, t_, q_, cp_, a_, v_, lambda_, delta_)
    T_t = (T_1_result + T_2_result) * 0.9
    return (T_t - T_n)/(T_1_result+T_2_result-T_n)
