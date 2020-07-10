# vs_circuit_solver.py
# номер версии не присвоен
# язык Python
#
# программа подбора значений R,C для вариантов электронной схемы
# исходя из моделирования подобной схемы в ngspice
# поставляется без всякой оптимизации, ибо имеет целью установление методики 
# расчета таких вещей и определения границ применимости этой методики
#
# автор В.Симонов, 30-июнь-2020
# vasily_simonov@mail.ru, github.com/vasily84
#
# license : это модуль в любом виде можно использовать в любых целях. 
# Ссылка на автора приветствуется, но не является обязательной 
#


import scipy.optimize as spo
import scipy.fft as scf
import math
import numpy as np
import matplotlib.pyplot as plt
from ctypes import c_double
import csv
import json
# внешние модули
import MySpice as spice
import libivcmp
import numpy


### SETTINGS ################################################################

# метод сравнения кривых тока и напряжения
# может быть : 'libivcmp','sko','sko_fft','power','power_fft'
#MISFIT_METHOD = 'sko'
MISFIT_METHOD = 'libivcmp'

# частота, Гц
INIT_F = 1e3
# амплитудное напряжение, Вольт, может изменится при загрузке внешнего файла данных
INIT_V = 2.5

# токоограничивающий резистор, Ом
INIT_Rcs = 4.7e-5

# SIGNAL/NOISE ratio
INIT_SNR = 135.0

# число циклов колебаний напряжения в записи
INIT_CYCLE = 3

# падение напряжения на диоде
# Диод считается полностью проводимым при напряжении больше чем DIODE_VOLTAGE,
# при меньшем полность закрыт. (Приближение)
DIODE_VOLTAGE = 0.6

# напряжение, при котором диоды считаем закрытыми
SMALL_VOLTAGE = 0.01

# "огромное сопротивление".
HUGE_R = 1e10 # 10 ГОм

# "большое сопротивление"
BIG_R = 1e8 # 100 МОм
 
# "мизерное сопротивление"
NULL_R = 1e-6 # 1 мкОм

# "малое сопротивление"
SMALL_R = 1e-4 # 100 мкОм

# "мизерная емкость"
NONE_C = 1e-15 # 0.001 пФ
 
# погрешность подбора кривых- критерий остановки. Подбор длится до тех пор, 
# пока функция сравнения не вернет значение CompareIvc()<=IVCMP_TOLERANCE     
IVCMP_TOLERANCE = 5e-2

# погрешность подбора номиналов в процентах. Номиналы емкостей считаются по 
# реактивному сопротивлению!. Подробности см. scipy.minimize(method='Powell')
VALUES_TOLERANCE = 1e-2

# число вычислений функции в процессе оптимизации. При малых значениях-
# минимально возможное число
MAXFEV = 100

# число точек в массивах тока и напряжения, может измениться при загрузке 
# внешнего файла данных
MAX_NUM_POINTS = 100


min_ivc = 1 
#############################################################################

# результат последнего моделирования в PySpice
analysis = None

# целевая кривая с током. Та, которую мы подбираем.
target_VCurrent = None
target_input_dummy = None

target_fileName = ''

# целевая кривая с током для сравнения в библиотеке libivcmp
target_IVCurve = None


# название временного файла схемы для запуска PySpice
circuit_SessionFileName = 'var1.cir'

# список значений для файла шаблона схемы. Число элементов - не меньше, чем 
# знаков {} в файле шаблона схемы
Xi_long = [0.,0.,0.,0., 0.,0.,0., 0.,0.,0.,0.]


# Маска оптимизируемых параметров - список булевого типа, например -
# Xi_long = [a, b, c, d]
# Xi_mask = [False,True,False,True] -> X_short = [b,d]
Xi_mask = [False,False,False,False, False,False,False, False,False,False,False]


#### ФУНКЦИИ ДЛЯ ШАБЛОНА, ЦЕЛЕВОЙ МОДЕЛИ И МАСКИ ПАРАМЕТРОВ ##################
def Xi_unroll(x_short):
    XL = Xi_long.copy()
    
    j = 0
    for i in range(0,len(Xi_mask)):
        if Xi_mask[i]:
            XL[i] += x_short[j] # 
            j += 1
                   
    return XL

def Xi_pack(Xi_):
    xi = []
       
    for i in range(0,len(Xi_mask)):
        if Xi_mask[i]:
            xi += [Xi_[i]]
                       
    return xi
    

# установить все известные номиналы
def set_circuit_nominals(nominals):
    global Xi_long  
    Xi_long = nominals.copy()
    
    
def reset_Xi_variable():
    for i in range(len(Xi_mask)):
        Xi_mask[i] = False
       
        
def set_Xi_variable(vlist):
    for v in vlist:
        if v=='R1': Xi_mask[0] = True
        if v=='C1': Xi_mask[1] = True
        if v=='_R_C1': Xi_mask[2] = True
        if v=='_R_D1': Xi_mask[3] = True
        
        if v=='R2': Xi_mask[4] = True
        if v=='C2': Xi_mask[5] = True
        if v=='_R_C2': Xi_mask[6] = True
        
        if v=='R3': Xi_mask[7] = True
        if v=='C3': Xi_mask[8] = True
        if v=='_R_C3': Xi_mask[9] = True
        if v=='_R_D3': Xi_mask[10] = True
    
        
# установить файл шаблона схемы. Считать его в набор строк.
def set_circuit_template(fileName):
    global circuit_TemplateStrings
    circuit_TemplateStrings = []
    templateF = open(fileName)
    
    for tStr in templateF:
        circuit_TemplateStrings += [tStr]
        #
    templateF.close()
   

# инициализировать целевую модель, промоделировав файл схемы
def init_target_by_circuitFile(fileName = circuit_SessionFileName):
    global target_VCurrent, target_input_dummy, target_IVCurve
    global circuit_SessionFileName
    
    var1 = circuit_SessionFileName
    circuit_SessionFileName = fileName
    process_circuitFile()
    circuit_SessionFileName = var1
    
    target_VCurrent = analysis.VCurrent
    target_input_dummy = analysis.input_dummy
    
    iv_curve = libivcmp.IvCurve()
    for i in range(MAX_NUM_POINTS):
        iv_curve.voltages[i] = c_double(analysis.input_dummy[i]) # Ток и напряжение были поменяны местами
        iv_curve.currents[i] = c_double(analysis.VCurrent[i])
    
    min_var_c = 0.01 * np.max(iv_curve.currents[:MAX_NUM_POINTS]) # value of noise for current
    min_var_v = 0.01 * np.max(iv_curve.voltages[:MAX_NUM_POINTS]) # value of noise for voltage
    libivcmp.SetMinVC(min_var_v, min_var_c) # Правильные значения фильтров для корректной работы
    
    target_IVCurve = iv_curve
              
           
# инициализировать целевую модель данными из csv файла, установить число точек на кривой MAX_NUM_POINTS
# определенными из файла 
def init_target_from_csvFile(fileName):
    global MAX_NUM_POINTS,INIT_V,NORMA_misfit
    global target_fileName
    global target_VCurrent,target_input_dummy,target_IVCurve
    #
    target_fileName = fileName
    
    with open(fileName,newline='') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')
        i = 0
        for row in csv_reader:
            # todo непонятно, почему читает пустые нечетные строки, выяснить
            if i==0:
                target_input_dummy = np.array(row,dtype=float)
            if i==2:
                target_VCurrent = np.array(row,dtype=float) 
            i += 1         
    
    MAX_NUM_POINTS = len(target_input_dummy)
    INIT_V = np.amax(target_input_dummy)

    iv_curve = libivcmp.IvCurve()
    for i in range(MAX_NUM_POINTS):
        iv_curve.voltages[i] = c_double(target_input_dummy[i]) # Ток и напряжение были поменяны местами
        iv_curve.currents[i] = c_double(target_VCurrent[i])
    
    min_var_c = 0.01 * np.max(iv_curve.currents[:MAX_NUM_POINTS]) # value of noise for current
    min_var_v = 0.01 * np.max(iv_curve.voltages[:MAX_NUM_POINTS]) # value of noise for voltage
    libivcmp.SetMinVC(min_var_v, min_var_c) # Правильные значения фильтров для корректной работы
    
    # libivcmp.set_MAX_NUM_POINTS(MAX_NUM_POINTS) # В некоторых случаях искажает работу модуля
    target_IVCurve = iv_curve
    
    try:
        NORMA_misfit=1./target_misfit_norma1()
    except ArithmeticError:
        NORMA_misfit=1.
        
    return


def Xi_to_RC(Xi):    
    RC = Xi.copy()
    
    RC[0] = np.abs(Xi[0])
    RC[1] = np.abs(R_to_C(Xi[1])) # C1
    RC[2] = np.abs(Xi[2])
    RC[3] = np.abs(Xi[3])
    
    RC[4] = np.abs(Xi[4])
    RC[5] = np.abs(R_to_C(Xi[5])) # C2
    RC[6] = np.abs(Xi[6])
    
    RC[7] = np.abs(Xi[7])
    RC[8] = np.abs(R_to_C(Xi[8])) # C3
    RC[9] = np.abs(Xi[9])
    RC[10] = np.abs(Xi[10])   
    return RC


# в наборе строк шаблона схемы сделать замену {} на значения 
# варьирования Xi_values, сохранить заданным с именем
def generate_circuitFile_by_values( Xi_values):
    newF = open(circuit_SessionFileName, 'w')
    rc_values = Xi_to_RC(Xi_values)
    
    newF.write('* cir file corresponding to the equivalent circuit.\n')
    # * Цепь 1
    if rc_values[0]<BIG_R: # цепь R1 присутствует
        if rc_values[2]>= BIG_R: # C1 присутствует 
            newF.write('R1 _net1 Input {:e}\n'.format(rc_values[0]))
            newF.write('C1 _net0 _net1 {:e}\n'.format(rc_values[1]))
        else: # С1 нет
            newF.write('R1 _net0 Input {:e}\n'.format(rc_values[0]))
        
        if rc_values[3]>= BIG_R: # D1 присутствует
            newF.write('D1 _net0 0 DMOD_D1 AREA=1.0 Temp=26.85\n')
        else: # вместо D1 перемычка
            newF.write('R_D1 0 _net0 {:e}\n'.format(rc_values[3]))
      
    # * Цепь 2
    if rc_values[4]<BIG_R:
        if rc_values[6]>= BIG_R: # C2 присутствует
            newF.write('R2 _net4 Input {:e}\n'.format(rc_values[4]))
            newF.write('C2 0 _net4 {:e}\n'.format(rc_values[5]))
        else: # вместо С2 перемычка, R2 сразу на землю
            newF.write('R2 0 Input {:e}\n'.format(rc_values[4]))
    
    # * Цепь 3
    if rc_values[7]<BIG_R:
        if rc_values[9]>=BIG_R: # C3 присутствует
            newF.write('R3 _net3 Input {:e}\n'.format(rc_values[7]))
            newF.write('C3 _net2 _net3 {:e}\n'.format(rc_values[8]))
        else: # С3 нет
            newF.write('R3 _net2 Input {:e}\n'.format(rc_values[7]))
            
        if rc_values[10]>=BIG_R: # D3 присутствует
            newF.write('D3 0 _net2 DMOD_D1 AREA=1.0 Temp=26.85\n')
        else: # вместо D3 перемычка
            newF.write('R_D3 0 _net2 {:e}\n'.format(rc_values[10]))
     
    # есть диоды, добавляем модель
    if (rc_values[10]>=BIG_R)or(rc_values[3]>= BIG_R):    
        newF.write('.MODEL DMOD_D1 D (Is=2.22e-10 N=1.65 Cj0=4e-12 M=0.333 Vj=0.7 Fc=0.5 Rs=0.0686 Tt=5.76e-09 Ikf=0 Kf=0 Af=1 Bv=75 Ibv=1e-06 Xti=3 Eg=1.11 Tcv=0 Trs=0 Ttt1=0 Ttt2=0 Tm1=0 Tm2=0 Tnom=26.85 )\n')
    
    newF.write('.END')     
    newF.close()

input_data = None
# промоделировать файл схемы
def process_circuitFile():
    global analysis,input_data
    #del(analysis) # необходимо
    
    if input_data is None:
        input_data = spice.Init_Data(INIT_F, INIT_V, INIT_Rcs,INIT_SNR )
        
    circuit = spice.LoadFile(circuit_SessionFileName)
    analysis = spice.CreateCVC1(circuit, input_data, MAX_NUM_POINTS, "input", INIT_CYCLE)   
    
        
# последний анализ перевести в форму, пригодную для сравнения в libivcmp
iv_curve = None
def analysis_to_IVCurve():
    global iv_curve
    
    if iv_curve is None:    
        iv_curve = libivcmp.IvCurve()
        
    for i in range(MAX_NUM_POINTS):
        iv_curve.voltages[i] = c_double(analysis.input_dummy[i])
        iv_curve.currents[i] = c_double(analysis.VCurrent[i])
        
    return iv_curve
    

def V_div_I(v,i):
    try:
        r = v/i
    except ArithmeticError:
        r = HUGE_R
    return r

    
    
# вывести на график результат моделирования
def analysis_plot(title='',pngName=''):
    plt.figure(1, (20, 10))
    plt.grid()
    
    # целевая ВАХ
    plt.plot(target_input_dummy, target_VCurrent,color='red')
    # ВАХ результат подбора
    plt.plot(analysis.input_dummy, analysis.VCurrent,color='blue')  
    
    s=''    
    if (not title==''):
       s = title
    elif not target_fileName=='':
       s = target_fileName
     

    s = s+', misfit='+ format(misfit_result,'0.5E')+', ivcmp='+format(ivcmp_result,'0.5E')
    
    plt.title(s)
    
    plt.xlabel('Напряжение [В]')
    plt.ylabel('Сила тока [А]')
    
    if(not pngName==''):
        plt.savefig(pngName)
        
    plt.show()
        

#### ФУНКЦИИ СРАВНЕНИЯ ВАХ ################################################### 
def C_to_R(c):
    r = 1/(2.*np.pi*INIT_F*c)
    return r


def R_to_C(r):
    c = 1/(2.*np.pi*INIT_F*r)
    return c

NORMA_misfit = 1.
def analysis_misfit():
    # суммирование без потери точности
    return NORMA_misfit*math.fsum(analysis_misfit_core()) 

def analysis_misfit_ivcmp():
    global min_ivc
    step_IVCurve = analysis_to_IVCurve()
    res = libivcmp.CompareIvc(target_IVCurve, step_IVCurve, MAX_NUM_POINTS)   
    if min_ivc > res:
        min_ivc = res
    #print(res, min_ivc)
    return res

# вычислить несовпадение последнего анализа и целевой функции.
def target_misfit_norma1():    
    curr_t = target_VCurrent
    volt_t = target_input_dummy
    
    if MISFIT_METHOD == 'power':
        r = (curr_t*volt_t)
        return math.fsum(r)
    
    if MISFIT_METHOD == 'power_fft':
        r = scf.rfft(curr_t*volt_t)       
        return math.fsum(r)
    
    if MISFIT_METHOD == 'sko':
        r = (curr_t)        
        r2 = np.abs(r)
        return math.fsum(r2)
    
    if MISFIT_METHOD == 'sko_fft':
        r = (curr_t)
        r2 = scf.rfft(np.abs(r))
        return math.fsum(r2)
    
    if MISFIT_METHOD == 'libivcmp':                
        return 1. # для libivcmp норма не поддерживается
    
    s = "unknown MISFIT_METHOD = '"+str(MISFIT_METHOD)+"'"
    raise RuntimeError(s)


# вычислить несовпадение последнего анализа и целевой функции.
def analysis_misfit_core():    
    curr_t = target_VCurrent
    curr_a = analysis.VCurrent
    volt_t = target_input_dummy
    volt_a = analysis.input_dummy
    
    if MISFIT_METHOD == 'power':
        r = (curr_t*volt_t-curr_a*volt_a)
        return r
    
    if MISFIT_METHOD == 'power_fft':
        r = scf.rfft(curr_t*volt_t-curr_a*volt_a)       
        return r
    
    if MISFIT_METHOD == 'sko':
        r = (curr_t-curr_a)        
        r2 = np.abs(r)
        return r2
    
    if MISFIT_METHOD == 'sko_fft':
        r = (curr_t-curr_a)
        r2 = np.abs(r)
        return scf.rfft(r2)
    
    if MISFIT_METHOD == 'libivcmp':                
        step_IVCurve = analysis_to_IVCurve()
        res = libivcmp.CompareIvc(target_IVCurve, step_IVCurve, MAX_NUM_POINTS)
        r = np.array([res])
        return r
    
    s = "unknown MISFIT_METHOD = '"+str(MISFIT_METHOD)+"'"
    raise RuntimeError(s)
    
    
       
#### ФУНКЦИИ РЕШАТЕЛЯ ########################################################
# вектор, обеспечивающий минимум оптимизируемой функции
Xi_result = [0.,0.,0.,0., 0.,0.,0., 0.,0.,0.,0.] 

# текущий найденный минимум оптимизируемой функции 
misfit_result = 0.
# результат сравнения найденного минимума по функцией CompareIvc()
ivcmp_result = 0.

# счетчик числа вызовов функции оптимизатором
FitterCount = 0
BestMisfitCount = 0
FITTER_SUCCESS = False

def calculate_misfit(Xi):
    generate_circuitFile_by_values(Xi)
    process_circuitFile()
    misfit = analysis_misfit()
    return misfit

# функция вызывается оптимизатором
def fitter_subroutine(Xargs):
    global Xi_result,misfit_result,FitterCount,BestMisfitCount,ivcmp_result,FITTER_SUCCESS
    FitterCount += 1
    xi = Xi_unroll(Xargs)
    misfit = calculate_misfit(xi)
    
    #print("fCount="+str(FitterCount)+', misfit='+str(Mscalar)+', Xargs='+str(Xargs))
    # первый запуск
    if FitterCount<=1:
        Xi_result = xi.copy()
        misfit_result = misfit
        ivcmp_result = analysis_misfit_ivcmp()
        BestMisfitCount = 0
        #print("fCount="+str(FitterCount)+', mCount='+str(BestMisfitCount)+', misfit='+str(Mscalar)+', Xargs='+str(Xargs))
            
    # лучший случай
    if misfit<misfit_result:
        Xi_result = xi.copy()
        misfit_result = misfit
        ivcmp_result = analysis_misfit_ivcmp()
        BestMisfitCount += 1
        #print("fCount="+str(FitterCount)+', mCount='+str(BestMisfitCount)+', misfit='+str(Mscalar)+', Xargs='+str(Xargs))
    
    # дополнительная проверка
    if ivcmp_result<=IVCMP_TOLERANCE: # достигли необходимой точности       
        FITTER_SUCCESS = True
        
    return misfit 

# 
def fitter_callback(Xk):
    global FITTER_SUCCESS
    if ivcmp_result<=IVCMP_TOLERANCE: # достигли необходимой точности       
        FITTER_SUCCESS = True
        return True
    
    return False    


# запустить автоподбор - сравнение по сумме отклонений точек         
def run_fitter(result_cir_file_name='',result_csv_file_name=''):                 
    Xargs = Xi_pack(Xi_long)
    
    for i in range(0,len(Xargs)):
        Xargs[i] = 0.
    
    resX = spo.minimize(fitter_subroutine,Xargs,method='Powell',callback=fitter_callback, options={'maxfev':MAXFEV,'xtol':VALUES_TOLERANCE})    
    #resX = spo.minimize(fitter_subroutine,Xargs,method='TNC',callback=fitter_callback, options={'maxfev':MAXFEV,'xtol':VALUES_TOLERANCE})    
              
    if(not result_csv_file_name==''):
        spice.SaveFile(analysis, result_csv_file_name)
    if(not result_cir_file_name==''):          
        generate_circuitFile_by_values(resX.x)
    
    return True
    
### элементарная схема ###
##############################################################################
def Sch_init():
    R = 100.
    C = 1e-6
    sch = {}
    sch['R1'] = R
    sch['C1'] = C # [1]
    sch['_R_C1'] = NULL_R 
    sch['_R_D1'] = HUGE_R
    
    sch['R2'] = R 
    sch['C2'] = C # [5]
    sch['_R_C2'] = NULL_R
    
    sch['R3'] = R 
    sch['C3'] = C #[8]
    sch['_R_C3'] = NULL_R 
    sch['_R_D3'] = HUGE_R         
    return sch


def Sch_get_Xi(sch):
    xi = []
    for k in sch:
        if(k=='C1')or(k=='C2')or(k=='C3'):
            xi+= [C_to_R(sch[k])]
        else:
            xi += [sch[k]]
    
    return xi


def Sch_load_from_Xi(sch,Xi):
    j = 0
    for k in sch:
        if(k=='C1')or(k=='C2')or(k=='C3'):
            sch[k] = R_to_C(Xi[j])
        else:
            sch[k] = Xi[j]
        j += 1

    
def get_Z1phase_approximation():
    # R1,C1 по максимуму тока и напряжения
    I_max = np.amax(target_VCurrent)
    max_cur_ind = np.argmax(target_VCurrent)
    max_volt_ind = np.argmax(target_input_dummy)
    
    phase_1 = 2.*np.pi*(max_cur_ind-max_volt_ind)/(MAX_NUM_POINTS)
    
    V_max = target_input_dummy[max_cur_ind]
    Z1 = V_div_I(V_max,I_max)
    return (Z1,phase_1)


def get_Z3phase_approximation():
    I_min = np.amin(target_VCurrent)
    
    # R3,C3 по минимуму тока и напряжения
    min_cur_ind = np.argmin(target_VCurrent)
    min_volt_ind = np.argmin(target_input_dummy)
    
    V_min = target_input_dummy[min_cur_ind]
    Z3 = V_div_I(V_min,I_min)
    
    phase_3 = 2.*np.pi*(min_cur_ind-min_volt_ind)/(MAX_NUM_POINTS)
    
    # R3 = Z3*np.cos(phase_3)
    # C3 = R_to_C(Z3*np.sin(phase_3))     
    return (Z3,phase_3)
    
    
def get_Z1phase_diode_approximation():
    # R1,C1 по максимуму тока и напряжения
    I_max = np.amax(target_VCurrent)
    max_cur_ind = np.argmax(target_VCurrent)
    max_volt_ind = np.argmax(target_input_dummy)
    
    phase_1 = 2.*np.pi*(max_cur_ind-max_volt_ind)/(MAX_NUM_POINTS)
    
    V_max = target_input_dummy[max_cur_ind]-DIODE_VOLTAGE
    Z1 = V_div_I(V_max,I_max)
    return (Z1,phase_1)


def get_Z3phase_diode_approximation():
    I_min = np.amin(target_VCurrent)
    
    # R3,C3 по минимуму тока и напряжения
    min_cur_ind = np.argmin(target_VCurrent)
    min_volt_ind = np.argmin(target_input_dummy)
    
    V_min = target_input_dummy[min_cur_ind]+DIODE_VOLTAGE
    Z3 = V_div_I(V_min,I_min)
    
    phase_3 = 2.*np.pi*(min_cur_ind-min_volt_ind)/(MAX_NUM_POINTS)
    
    # R3 = Z3*np.cos(phase_3)
    # C3 = R_to_C(Z3*np.sin(phase_3))     
    return (Z3,phase_3)


CODE2_COUNT = 4
# ses - сессия варьирования, которую необходимо проинициализировать
# swcode - числовой код,от 0 до 255 включительно, задает положения переключателей
# code2 - дополнительный код, для каждого варианта swcode передавать code2=0,1,2, ...
# до тех пор, пока функция не вернет False
def Session_init_by_approximation(ses,swcode,code2):    
    R1=False if(swcode & 1) else True    
    R2=False if(swcode & 2) else True
    R3=False if(swcode & 4) else True    
    # C1=False if(swcode & 8) else True         
    # D1=False if(swcode & 16) else True     
    # C3=False if(swcode & 32) else True   
    # D3=False if(swcode & 64) else True             
    # C2=False if(swcode & 128) else True   
    
    sch = ses['start_sch']           
    Sch_init_by_approximation(sch)
    Session_set_switchers(ses,swcode)
    Session_run1(ses)
    return False

    if code2==0:
        Session_run1(ses)
        return True # вызвать снова

    if R2 and (not R1)and(not R3):
        (r2,phase2) = get_Z1phase_approximation()
        sch['R2']=r2*np.cos(np.pi*code2/CODE2_COUNT)
        #sch['C2']=R_to_C(r2*np.sin(np.pi*code2/CODE2_COUNT))
        ses['start_sch'] = sch     
        Session_run1(ses)
        #if code2>=CODE2_COUNT: return False
        return False
    
    if R2 and R1 and (not R3):
        (r2,phase2) = get_Z3phase_approximation()
        sch['R2']=r2*np.cos(phase2)
        #sch['C2']=R_to_C(r2*np.sin(phase2))
        
        (r1,phase1) = get_Z1phase_diode_approximation()
        sch['R1']=r1*np.cos(phase1)
        #sch['C1']=R_to_C(r1*np.sin(phase1))
        ses['start_sch'] = sch
        Session_run1(ses)
        return False
    
    if R2 and R3 and (not R1):
        (r2,phase2) = get_Z1phase_approximation()
        sch['R2']=r2*np.cos(phase2)
        #sch['C2']=R_to_C(r2*np.sin(phase2))
        
        (r3,phase3) = get_Z3phase_diode_approximation()
        sch['R3']=r3*np.cos(phase3)
        #sch['C3']=R_to_C(r3*np.sin(phase3))
        ses['start_sch'] = sch
        Session_run1(ses)
        return False
    
    if R1 and R3 and (not R2):
        (r1,phase1) = get_Z1phase_diode_approximation()
        sch['R1']=r1*np.cos(phase1)
        #sch['C1']=R_to_C(r1*np.sin(phase1))
        
        (r3,phase3) = get_Z3phase_diode_approximation()
        sch['R3']=r3*np.cos(phase3)
        #sch['C3']=R_to_C(r3*np.sin(phase3))
        ses['start_sch'] = sch
        Session_run1(ses)
        return False
    
    if R1 and (not R2)and(not R3):
        (r1,phase1) = get_Z1phase_diode_approximation()
        sch['R1']=r1*np.cos(phase1)
        #sch['C1']=R_to_C(r1*np.sin(phase1))
        ses['start_sch'] = sch
        Session_run1(ses)
        return False
    
    if R3 and (not R1)and(not R2):
        (r3,phase3) = get_Z3phase_diode_approximation()
        sch['R3']=r3*np.cos(phase3)
        #sch['C3']=R_to_C(r3*np.sin(phase3))
        ses['start_sch'] = sch
        Session_run1(ses)
        return False
    
    Session_run1(ses)
    return False


def Sch_init_by_approximation(sch):
    # R1,C1 по максимуму тока и напряжения
    I_max = np.amax(target_VCurrent)
    max_cur_ind = np.argmax(target_VCurrent)
    max_volt_ind = np.argmax(target_input_dummy)
    
    phase_1 = 2.*np.pi*(max_cur_ind-max_volt_ind)/(MAX_NUM_POINTS)
    
    V_max = target_input_dummy[max_cur_ind]-DIODE_VOLTAGE
    Z1 = V_div_I(V_max,I_max)
    #Z1_ = V_div_I(V_max+DIODE_VOLTAGE,I_max)
    R1 = Z1*np.cos(phase_1)
    C1 = R_to_C(Z1*np.sin(phase_1))
    
    sch['R1'] = np.abs(R1)
    if np.abs(phase_1*180./np.pi)>1.:
        sch['_R_C1'] = HUGE_R
        sch['C1'] = C1
    
    # print('R1='+str(R1)+', C1='+str(C1)+' Z1='+str(Z1)+' ,Z1_='+str(Z1_))  
    # print('phase_1='+str(phase_1*180./np.pi))
    I_min = np.amin(target_VCurrent)
    
    # R3,C3 по минимуму тока и напряжения
    min_cur_ind = np.argmin(target_VCurrent)
    min_volt_ind = np.argmin(target_input_dummy)
    
    V_min = target_input_dummy[min_cur_ind]+DIODE_VOLTAGE
    Z3 = V_div_I(V_min,I_min)
    #Z3_ = V_div_I(V_min-DIODE_VOLTAGE,I_min)
    
    phase_3 = 2.*np.pi*(min_cur_ind-min_volt_ind)/(MAX_NUM_POINTS)
    
    R3 = Z3*np.cos(phase_3)
    C3 = R_to_C(Z3*np.sin(phase_3))
    # print('R3='+str(R3)+', C3='+str(C3)+', Z3='+str(Z3)+', Z3_='+str(Z3_))
    # print('phase_3='+str(phase_3*180./np.pi))
    
    sch['R3'] = np.abs(R3)
    if np.abs(phase_3*180./np.pi)>1.:
        sch['_R_C3'] = HUGE_R
        sch['C3'] = C3
        
    R2 = R1+R3  # переделать на окрестности малого сигнала
    C2 = C1+C3  # переделать на сдвиг фаз в окрестности ноля
    
    sch['R2'] = R2
    if C_to_R(C2)>0.01*R2:
        sch['_R_C3'] = HUGE_R
        sch['C2'] = C2
        
        
def Sch_init_by_approximation_old(sch):
    # R1,C1 по максимуму тока и напряжения
    I_max = np.amax(target_VCurrent)
    max_cur_ind = np.argmax(target_VCurrent)
    max_volt_ind = np.argmax(target_input_dummy)
    
    phase_1 = 2.*np.pi*(max_cur_ind-max_volt_ind)/(MAX_NUM_POINTS)
    
    V_max = target_input_dummy[max_cur_ind]-DIODE_VOLTAGE
    Z1 = V_div_I(V_max,I_max)
    #Z1_ = V_div_I(V_max+DIODE_VOLTAGE,I_max)
    R1 = Z1*np.cos(phase_1)
    C1 = R_to_C(Z1*np.sin(phase_1))
    
    sch['R1'] = np.abs(R1)
    if np.abs(phase_1*180./np.pi)>1.:
        sch['_R_C1'] = HUGE_R
        sch['C1'] = C1
    
    # print('R1='+str(R1)+', C1='+str(C1)+' Z1='+str(Z1)+' ,Z1_='+str(Z1_))  
    # print('phase_1='+str(phase_1*180./np.pi))
    I_min = np.amin(target_VCurrent)
    
    # R3,C3 по минимуму тока и напряжения
    min_cur_ind = np.argmin(target_VCurrent)
    min_volt_ind = np.argmin(target_input_dummy)
    
    V_min = target_input_dummy[min_cur_ind]+DIODE_VOLTAGE
    Z3 = V_div_I(V_min,I_min)
    #Z3_ = V_div_I(V_min-DIODE_VOLTAGE,I_min)
    
    phase_3 = 2.*np.pi*(min_cur_ind-min_volt_ind)/(MAX_NUM_POINTS)
    
    R3 = Z3*np.cos(phase_3)
    C3 = R_to_C(Z3*np.sin(phase_3))
    # print('R3='+str(R3)+', C3='+str(C3)+', Z3='+str(Z3)+', Z3_='+str(Z3_))
    # print('phase_3='+str(phase_3*180./np.pi))
    
    sch['R3'] = np.abs(R3)
    if np.abs(phase_3*180./np.pi)>1.:
        sch['_R_C3'] = HUGE_R
        sch['C3'] = C3
        
    R2 = R1+R3  # переделать на окрестности малого сигнала
    C2 = C1+C3  # переделать на сдвиг фаз в окрестности ноля
    
    sch['R2'] = R2
    if C_to_R(C2)>0.01*R2:
        sch['_R_C3'] = HUGE_R
        sch['C2'] = C2
        
def Sch_saveToFile(sch,fileName):
    try:
        newF = open(fileName, 'w')
        json.dump(sch,newF)
    finally:           
        newF.close()
    
def init_target_by_Sch(sch):
    global NORMA_misfit
    generate_circuitFile_by_values(Sch_get_Xi(sch))
    init_target_by_circuitFile()  
    
    try:
        NORMA_misfit=1./target_misfit_norma1()
    except ArithmeticError:
        NORMA_misfit=1.
            
        
############################################################################
## площадь кривой X(i),Y(i), в пределах индексов
def curve_SquareXY(X,Y,index_from,index_to):
    S = 0.
    for i in range(index_from,index_to):
        dx = X[i+1]-X[i]
        y = (Y[i+1]+Y[i])/2.
        S += dx*y
    
    return S
        

def curve_INDEX(voltage,current):
    fsm_region = 0
    i0=0;i1=0;i2=0;i3=0;i4=0;i5=0;i6=0; i7 = len(voltage)-1
    
    L2 = voltage**2+current**2 # квадрат длины
    max_L2 = 0
    
    for i in range(0,len(voltage)):
        if fsm_region==0:
            if voltage[i]<(-DIODE_VOLTAGE-SMALL_VOLTAGE):
                i1 = i
                i2 = i
                max_L2 = L2[i]
                fsm_region = 1
                continue
            
        if fsm_region==1:
            if L2[i]>max_L2:
                max_L2 = L2[i]
                i2 = i
                
            if voltage[i]>(-DIODE_VOLTAGE+SMALL_VOLTAGE):
                i3 = i
                fsm_region = 2
                continue
            
        if fsm_region==2:
            if voltage[i]>(DIODE_VOLTAGE+SMALL_VOLTAGE):
                i4 = i
                i5 = i
                max_L2 = L2[i]
                fsm_region = 3
                continue
            
        if fsm_region==3:
            if L2[i]>max_L2:
                max_L2 = L2[i]
                i5 = i
                
            if voltage[i]<(DIODE_VOLTAGE-SMALL_VOLTAGE):
                i6 = i
                fsm_region = 7
                continue
    ## end_for
    
    S1 = curve_SquareXY(voltage,current,i0,i1)
    
    S2 = curve_SquareXY(voltage,current,i1,i2)
    S3 = curve_SquareXY(voltage,current,i2,i3)
    
    S4 = curve_SquareXY(voltage,current,i3,i4)
    
    S5 = curve_SquareXY(voltage,current,i4,i5)
    S6 = curve_SquareXY(voltage,current,i5,i6)
    
    S7 = curve_SquareXY(voltage,current,i6,i7)
    
    
    print('i0='+str(i0)+', i1='+str(i1)+', i2='+str(i2)+', i3='+str(i3)+', i4='+str(i4)+', i5='+str(i5)+' ,i6='+str(i6)+', i7='+str(i7))
    print('S1='+str(S1)+', S2='+str(S2)+', S3='+str(S3)+', S4='+str(S4)+', S5='+str(S5)+', S6='+str(S6)+' ,S7='+str(S7))
    
    return (i0,i1,i2,i3,i4,i5,i6,i7)      

# вычислить фазовый сдвиг для уровней сигнала меньше напряжения диода на 
# отрицательной полуволне
def calc_phase_shift1():
    pass

def calc_phase_shift2():
    pass

def calc_phase_shift3():
    pass
                
    
############################################################################# 
def Session_create(start_sch):
    s = {}
    s['start_sch'] = start_sch
    return s

# выполнить схему один раз
def Session_run1(session):
    try:
        sch = session['result_sch']
    except KeyError:
        sch = session['start_sch']
    
    #print('\nSession_run1: sch='+str(sch))
    xi = Sch_get_Xi(sch)
    set_circuit_nominals(xi)
    session['misfit'] = calculate_misfit(xi)
    #print('Session_run1: misfit='+str(session['misfit']))
        
    
# запустить подбор для сессии
def Session_run_fitter(session):
    global FitterCount
    FitterCount = 0
    try:
        sch = session['result_sch']
    except KeyError:
        sch = session['start_sch']
    else:
        session['start_sch']=sch
        
    var_list = session['Xi_variable']
    set_circuit_nominals(Sch_get_Xi(sch))
    set_Xi_variable(var_list)
    run_fitter()
    sch2 = Sch_init()
    Sch_load_from_Xi(sch2, Xi_result)
    session['result_sch'] = sch2
    
    session['misfit'] = misfit_result
    session['fCount'] = FitterCount
    session['mCount'] = BestMisfitCount
   

# проверить, имеет ли смысл такая установка переключателей в схеме
def is_valid_switchers(swcode):  
    if swcode & (1+2+3): # все ветви заглушены
        return False
    
    # все разыгрывание по заглушенной первой ветке
    if (swcode==1)or(swcode==1+8)or(swcode==1+16)or(swcode==1+8+16):
        return False

    # все разыгрывание по заглушенной второй ветке
    if (swcode==2)or(swcode==2+128):
        return False
    
    # все разыгрывание по заглушенной третьей ветке
    if (swcode==3)or(swcode==3+32)or(swcode==3+64)or(swcode==3+32+64):
        return False
    
    return True


# установить переключатели для схемы.
def Session_set_switchers(session, swcode):
    sch = session['start_sch']
    var_list = []
    
    if swcode & 1: # ветка 1  
        #print('dnp branch1')
        sch['R1'] = HUGE_R
    else:
        var_list +=['R1']
        
    if swcode & 2: # ветка 2 
        #print('dnp branch2') 
        sch['R2'] = HUGE_R
    else:
        var_list += ['R2']
    
    if swcode & 4: # ветка 3  
        #print('dnp branch3')
        sch['R3'] = HUGE_R
    else:
        var_list += ['R3']
        
    if swcode & 8: # C1  
        #print('dnp C1')
        sch['_R_C1'] = NULL_R
    else:
        sch['_R_C1'] = HUGE_R
        var_list += ['C1']
        
    if swcode & 16: # D1 
        #print('dnp D1') 
        sch['_R_D1'] = NULL_R
    else:
        sch['_R_D1'] = HUGE_R
        
    
    if swcode & 32: # C3 
        #print('dnp C3')
        sch['_R_C3'] = NULL_R
    else:
        sch['_R_C3'] = HUGE_R
        var_list += ['C3']
        
    if swcode & 64: # D3  
        #print('dnp D3')
        sch['_R_D3'] = NULL_R
    else:
        sch['_R_D3'] = HUGE_R
        
    if swcode & 128: # C2
        #print('dnp C2')
        sch['_R_C2'] = NULL_R
    else:
        sch['_R_C2'] = HUGE_R
        var_list += ['C2']
    
    session['Xi_variable'] = var_list  
      
    
def Session_processAll(fileName='result.txt'):
    global FITTER_SUCCESS,VALUES_TOLERANCE,MAXFEV
    FITTER_SUCCESS = False
    ses_list = []
    best_ses = None
    best_misfit = 2 # заведомо большое число
    
    ses_count = 0
    ## создаем список сессий для старта
    for swcode in range(255):
        if not is_valid_switchers(swcode): continue
        
        code2 = 0    
        next_code2 = True
        
        while next_code2:
            sch0 = Sch_init()
            ses = Session_create(sch0)
            next_code2 = Session_init_by_approximation(ses,swcode,code2)
            code2 += 1
            ses_list += [ses]
            ses_count += 1
            
            if ses['misfit']<IVCMP_TOLERANCE: # условие останова удовлетворено
                best_ses = ses
                best_misfit = best_ses['misfit']
                print(ses['start_sch'])           
                analysis_plot('FITTER SUCCESS')
                print('FITTER_SUCCESS!!\nmisfit = '+str(best_misfit))
                Sch_saveToFile(ses['start_sch'], fileName)
                print('good case!!')
                return
            
    ## end_for   
    #print('ses_count = '+str(ses_count))

    ## сортируем сессии, чтобы начать подбор с наиболее подходящих
    ses_list = sorted(ses_list,key=lambda s:s['misfit'])  
    best_ses = ses_list[0]
    best_misfit = best_ses['misfit']
        
    ## запускаем автоподбор, пока не будут удовлетворены условия останова
    for ses in ses_list:
        Session_run_fitter(ses)
        if(ses['misfit']<best_misfit):
            best_misfit = ses['misfit']
            best_ses = ses
            print('misfit = '+str(best_misfit))
            #print(ses['result_sch'])           
            #analysis_plot()
            if FITTER_SUCCESS:
                print(ses['result_sch'])           
                analysis_plot('FITTER SUCCESS')
                print('FITTER_SUCCESS!!\nmisfit = '+str(best_misfit))
                Sch_saveToFile(ses['result_sch'], fileName)
                return
    ## end_for
    
    # подбор завершился неудачно, выводим что есть
    print('FITTER routine unsuccessfull\nmisfit = '+str(best_misfit))   
    Sch_saveToFile(best_ses['result_sch'], fileName)
    Session_run1(best_ses)
    analysis_plot('FITTER routine unsuccessfull')


#############################################################################
        
def test1():
    sch = Sch_init()
    sch['R1'] = 1e3
    sch['R2'] = 1e2
    init_target_by_Sch(sch)
    print('test1()')
    Session_processAll('test1.txt')
    

def test2():
    sch = Sch_init()
    sch['R1'] = 1e2
    sch['C1'] = 1e-5
    sch['_R_C1'] = HUGE_R
    sch['R3'] = 1e3
    init_target_by_Sch(sch)
    print('test2()')
    Session_processAll('test2.txt')

def test3():
    sch = Sch_init()
    sch['R2'] = 1e2
    sch['_R_C2'] = HUGE_R
    sch['C2'] = 1e-7
    init_target_by_Sch(sch) 
    print('test3()')
    Session_processAll('test3.txt')
    

def test4(): 
    sch = Sch_init()
    sch['R2'] = 1e2
    sch['_R_C2'] = HUGE_R
    sch['C2'] = 1e-7
    sch['R3'] = 1e3
    init_target_by_Sch(sch) 
    print('test4()')
    Session_processAll('test4.txt')
    
def test5(): 
    sch = Sch_init()
    sch['R1'] = 1e2
    sch['C1'] = NONE_C
    sch['_R_C1'] = NULL_R
    sch['_R_D1'] = NULL_R
    
    sch['R2'] = NULL_R
    sch['_R_C2'] = HUGE_R
    sch['C2'] = 1e-7
    init_target_by_Sch(sch) 
    print('test5()')
    Session_processAll('test5.txt')
    
def test_data(csv_data,fileName='result.txt'):
    print('\n')
    print(csv_data)
    init_target_from_csvFile(csv_data)
    plt.plot(target_VCurrent)
    plt.show()
    plt.plot(target_input_dummy)
    plt.show()
    Session_processAll(fileName)


def test_circuit(circuitFile,resultFile = 'result.txt'):
    print('\n')
    print(circuitFile)
    init_target_by_circuitFile(circuitFile)
    #
    # plt.plot(target_input_dummy)
    # plt.title('voltage')
    # plt.show()
    # plt.plot(target_VCurrent)
    # plt.title('current')
    # plt.show()
    # plt.plot(target_VCurrent*target_input_dummy)
    # plt.title('power')
    # plt.show()
    #curve_INDEX(target_input_dummy,target_VCurrent)
    Session_processAll(resultFile)
 
def main(): 
    #libivcmp.set_MAX_NUM_POINTS(100)
    #init_target_from_csvFile('пример1.csv')
    # test1()
    # test2()
    # test3()
    # test4()
    # test5()
    # return
    # test_data('test_nores.csv')
    # test_data('test_rc.csv')
    # test_data('test_rcd.csv')
    # test_data('test_rcd_no.csv')
    # test_data('test_rcd_rcd.csv')
    # test_data('test_rd_ii_c.csv')
    # return 
    test_circuit('NoRes.cir')
    test_circuit('R-C.cir')
    test_circuit('R-C-D.cir')
    test_circuit('R-C-D_No.cir')
    test_circuit('R-D_ii_C.cir')
    test_circuit('R1C1D1_ii_R2C2D2.cir')
    
    
if __name__=='__main__':
    main()
      
##############################################################################
  
