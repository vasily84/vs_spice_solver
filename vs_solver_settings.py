# vs_solver_settings.py
# номер версии не присвоен
# язык Python
#
# составная часть проекта vs_circuit_solver
# глобальные переменные для управления решателям
#
# автор В.Симонов, 18-мая-2020
# vasily_simonov@mail.ru, github.com/vasily84
#
# license : это модуль в любом виде можно использовать в любых целях. 
# Ссылка на автора приветствуется, но не является обязательной 
#


# использовать внешнюю библиотеку libivcmp для сравнения кривых?
# !! установить True или False 
USE_LIBIVCMP = False

# метод оптимизации функции подбора параметров R,C,
# варианты для функции scipy.optimize.minimum() 
# !! Раскомментировать необходимый FITTER_METHOD
FITTER_METHOD = 'Powell' # это метод работает лучше всего
#FITTER_METHOD = 'Nelder-Mead' # тоже рабочий
#FITTER_METHOD = None # не рабочий - происходит потеря точности 


# частота, Гц
INIT_F = 1000.
# амплитудное напряжение, В
INIT_V = 3.3
# токоограничивающий резистор, Ом
INIT_Rcs = .47

# SIGNAL/NOISE ratio
INIT_SNR = 40.0

# число циклов колебаний напряжения в записи
INIT_CYCLE = 1

# падение напряжения на диоде
# Диод считается полностью проводимым при напряжении больше чем DIODE_VOLTAGE,
# при меньшем полность закрыт. (Приближение)
DIODE_VOLTAGE = 0.7

# напряжение, при котором диоды считаем закрытыми
SMALL_VOLTAGE = 0.1



# "огромное сопротивление".
HUGE_R = 1e10 # 

# "большое сопротивление"
LARGE_R = 1e7 # 100 МегаОм
# "мизерное сопротивление"
NULL_R = 1e-6 # todo

# "мизерная емкость"
NONE_C = 1e-15 # 0.001 пФ
 
# погрешность подбора кривых
TOLERANCE = 1e-3
# число точек в массивах тока и напряжения

MAX_NUM_POINTS = 300


_TOLERANCE = TOLERANCE
_MAX_NUM_POINTS = MAX_NUM_POINTS 
_INIT_F = INIT_F
_INIT_V = INIT_V
_INIT_Rcs = INIT_Rcs
_INIT_SNR = INIT_SNR


def restore_INIT_param():
    global TOLERANCE,MAX_NUM_POINTS,INIT_F,INIT_V,INIT_Rcs,INIT_SNR
    TOLERANCE = _TOLERANCE
    MAX_NUM_POINTS = _MAX_NUM_POINTS 
    INIT_F = _INIT_F
    INIT_V = _INIT_V
    INIT_Rcs = _INIT_Rcs
    INIT_SNR = _INIT_SNR
