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

# метод сравнения кривых тока и напряжения

USE_LIBIVCMP = False
USE_FFT_FILTRATION = True
FFT_HIGH = 35 

#MISFIT_METHOD = 'libivcmp' # использовать внешнюю библиотеку libivcmp
MISFIT_METHOD = 'sko'
#MISFIT_METHOD = 'sko_fft'
#MISFIT_METHOD = 'power'
#MISFIT_METHOD = 'power_fft' 

MISFIT_KIND = 'minimize' # сравнение по сумме несовпадений в точках
#MISFIT_KIND = 'least_square' # сравнение по множеству точек


CIR_TEMPLATE = 'general.cir_t'

# метод оптимизации функции подбора параметров R,C,
# варианты для функции scipy.optimize.minimum() 
# !! Раскомментировать необходимый FITTER_METHOD
FITTER_METHOD = 'Powell' # это метод работает лучше всего
#FITTER_METHOD = 'Nelder-Mead' # тоже рабочий
#FITTER_METHOD = 'SLSQP'  


# частота, Гц
INIT_F = 1e3
# амплитудное напряжение, В
INIT_V = 2.5
# токоограничивающий резистор, Ом
INIT_Rcs = 0.47

# SIGNAL/NOISE ratio
INIT_SNR = 35.0

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
MAX_NUM_POINTS = 500



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
