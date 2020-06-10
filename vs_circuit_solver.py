# vs_circuit_solver.py
# номер версии не присвоен
# язык Python
#
# программа подбора значений R,C для вариантов электронной схемы
# исходя из моделирования подобной схемы в ngspice
# поставляется без всякой оптимизации, ибо имеет целью установление методики 
# расчета таких вещей и определения границ применимости этой методики
#
# автор В.Симонов, 18-мая-2020
# vasily_simonov@mail.ru, github.com/vasily84
#
# license : это модуль в любом виде можно использовать в любых целях. 
# Ссылка на автора приветствуется, но не является обязательной 
#

import scipy.optimize as spo
import scipy.fft as scf
import numpy as np
import matplotlib.pyplot as plt
from ctypes import c_double
import csv


# управляющие глобальные переменные
import vs_solver_settings as G
# модуль Михаила Лукьянова, github.com/LukyanovM/MySpice 
import MySpice as spice
# сравнение кривых - закомментировать, если этой библиотеки нет
import libivcmp
import random


if G.USE_LIBIVCMP:
    G.TOLERANCE = 1e-15
    G.MAX_NUM_POINTS = libivcmp.MAX_NUM_POINTS


# результат последненго моделирования в PySpice
analysis = None

# целевая кривая с током. Та, которую мы подбираем.
target_VCurrent = None
target_input_dummy = None


# целевая кривая с током для сравнения в библиотеке libivcmp
target_IVCurve = None


# название временного файла схемы для запуска PySpice
circuit_SessionFileName = 'var1.cir'

# строки - шаблон файла схемы PySpice. 
# Шаблон файла схемы - файл с расширением *.cir_t, аналогичный файлу описания 
# схемы ngspice *.cir, где вместо конкретного значения емкости или сопротивлния
# компонента стоит символ замены '{}'
circuit_TempateStrings = None

# список значений для файла шаблона схемы. Число элементов - не меньше, чем 
# знаков {} в файле шаблона схемы
Xi_long = [0.,0.,0.,0., 0.,0.,0., 0.,0.,0.,0.]

# погрешность отпредения независимой переменной
#Xi_step = [1.,1.,1.,1., 1.,1.,1., 1.,1.,1.,1.]

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
    global circuit_TempateStrings
    circuit_TempateStrings = []
    templateF = open(fileName)
    
    for tStr in templateF:
        circuit_TempateStrings += [tStr]
        #
    templateF.close()
   
    
# инициализировать целевую модель, промоделировав файл схемы
def init_target_by_circuitFile(fileName = circuit_SessionFileName):
    global target_VCurrent, target_input_dummy, target_IVCurve,target_Q
    process_circuitFile()
    target_VCurrent = analysis.VCurrent
    target_input_dummy = analysis.input_dummy
    
    if G.USE_LIBIVCMP:
        target_IVCurve = analysis_to_IVCurve()
              
           
# инициализировать целевую модель данными из csv файла, с заданной частотой,
# амплитудой, и токоограничивающим резистором
def init_target_from_csvFile(fileName, F=G.INIT_F, V=G.INIT_V, Rcs=G.INIT_Rcs):
    global target_VCurrent,target_input_dummy,target_IVCurve
    #
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
            
    
    # 
    G.INIT_F = F
    G.INIT_V = V
    G.INIT_Rcs = Rcs
    
    if G.USE_LIBIVCMP:
        iv_curve = libivcmp.IvCurve()
        for i in range(G.MAX_NUM_POINTS):
            iv_curve.voltages[i] = c_double(target_VCurrent[i])
            iv_curve.currents[i] = c_double(target_input_dummy[i])
        
        libivcmp.SetMinVC(0, 0)
        target_IVCurve = iv_curve
    
    return


def Xi_to_RC(Xi):    
    RC = Xi.copy()
    RC[1] = R_to_C(Xi[1]) # C1
    RC[5] = R_to_C(Xi[5]) # C2
    RC[8] = R_to_C(Xi[8]) # C3
    
    return RC

    
# в наборе строк шаблона схемы сделать замену {} на значения 
# варьирования Xi_values, сохранить заданным с именем
def generate_circuitFile_by_values( Xi_values):
    fileCircuit = circuit_SessionFileName
    #print('\ngenerate_circuitFile_by_values()\n'+str(Xi_values)) 
    newF = open(fileCircuit, 'w')
    rc_values = Xi_to_RC(Xi_values)
    i = 0
    
    for tStr in circuit_TempateStrings:
        cStr = tStr
        if tStr.find('{}')>=0:
            cStr = tStr.format(str(np.abs(rc_values[i]))) # ставим абсолютные значения номиналов
            i +=1
            
        newF.write(cStr)
    newF.close()

        
# промоделировать файл схемы
def process_circuitFile(csvName=''):
    global analysis
    del(analysis) # необходимо
    fileName=circuit_SessionFileName
    
    circuit = spice.LoadFile(fileName)
    input_data = spice.Init_Data(G.INIT_F, G.INIT_V, G.INIT_Rcs,G.INIT_SNR )
    analysis = spice.CreateCVC1(circuit, input_data, G.MAX_NUM_POINTS, "input", G.INIT_CYCLE)   
    
    if(not csvName==''):
        spice.SaveFile(analysis, csvName)
 
    
def my_fft_filtration(data, high=G.FFT_HIGH):
    if G.USE_FFT_FILTRATION:
        f_data = scf.rfft(data)
        f_data[high:] = 0.
        i_data = scf.irfft(f_data)
        return i_data
    return data

# последний анализ перевести в форму, пригодную для сравнения в libivcmp
def analysis_to_IVCurve():
    iv_curve = libivcmp.IvCurve()
    for i in range(G.MAX_NUM_POINTS):
        iv_curve.voltages[i] = c_double(analysis.VCurrent[i])
        iv_curve.currents[i] = c_double(analysis.input_dummy[i])
        
    libivcmp.SetMinVC(0, 0)
    return iv_curve
    

def V_div_I(v,i):
    try:
        r = v/i
    except ArithmeticError:
        r = G.HUGE_R
    return r


def create_stat_series():
    u = {}
    u['X_summ'] = 0.
    u['X2_summ'] = 0.
    u['count'] = 0
    return u


def add_to_stat_series(u,Value):
    u['X_summ'] += Value
    u['X2_summ'] += Value*Value
    u['count'] +=1
    
def my_R_plot():
    Iarr = my_fft_filtration(target_VCurrent)
    Varr = my_fft_filtration(target_input_dummy)
    Rarr = Varr/Iarr
    plt.plot(Varr,Rarr)
    plt.show()
    
def my_Sigma_plot():
    Iarr = my_fft_filtration(target_VCurrent)
    Varr = my_fft_filtration(target_input_dummy)
    Sigma_arr = Iarr/Varr
    #plt.plot(Varr,Sigma_arr)
    plt.plot(Sigma_arr)
    plt.show()
    
def my_Power_plot():
    Iarr = my_fft_filtration(target_VCurrent)
    Varr = my_fft_filtration(target_input_dummy)
    Power_arr = Iarr*Varr
    #plt.plot(Varr,Sigma_arr)
    plt.plot(Power_arr)
    plt.show()
    
def my_PowerPrime_plot():
    Iarr = my_fft_filtration(target_VCurrent)
    Varr = my_fft_filtration(target_input_dummy)
    Power_arr = Iarr*Varr
    #plt.plot(Varr,Sigma_arr)
    prime_arr = Power_arr
    for i in range(0,len(Power_arr)-1):
        prime_arr[i] = Power_arr[i+1]-Power_arr[i]
    plt.plot(prime_arr)
    
    plt.show()
    
    
    



# вывести на график результат моделирования
def analysis_plot(title='',pngName=''):
    figure1 = plt.figure(1, (20, 10))
    plt.grid()
    
    # целевая ВАХ
    plt.plot(target_input_dummy, target_VCurrent,color='red')
    # ВАХ результат подбора
    plt.plot(analysis.input_dummy, analysis.VCurrent,color='blue')  
           
    
    if (not title==''):
       plt.title(title)       
    plt.xlabel('Напряжение [В]')
    plt.ylabel('Сила тока [А]')
    if(not pngName==''):
        plt.savefig(pngName)
        
    plt.show()
    
    
# вывести на график спектр сигнала
def analysis_plotFFT():
    figure = plt.figure(1,(20,10))
    plt.grid()
    
    S1 = scf.rfft(target_VCurrent)
    plt.plot(np.abs(S1)[:30],color='red')
    
    S2 = scf.rfft(analysis.VCurrent)
    plt.plot(np.abs(S2)[:30],color='blue')
    
    S = S2-S1
    plt.plot(np.angle(S)[:30],color='yellow')
    
    plt.show()
    

#### ФУНКЦИИ СРАВНЕНИЯ ВАХ ################################################### 
def C_to_R(c):
    r = 1/(2.*np.pi*G.INIT_F*c)
    return r


def R_to_C(r):
    c = 1/(2.*np.pi*G.INIT_F*r)
    return c


def analysis_misfit_scalar():
    return np.sum(analysis_misfit_core())


# вычислить несовпадение последнего анализа и целевой функции.
def analysis_misfit_core():    
    curr_t = my_fft_filtration(target_VCurrent)
    curr_a = my_fft_filtration(analysis.VCurrent)
    volt_t = my_fft_filtration(target_input_dummy)
    volt_a = my_fft_filtration(analysis.input_dummy)
    
    if G.MISFIT_METHOD == 'power':
        r = (curr_t*volt_t-curr_a*volt_a)
        return r
    
    if G.MISFIT_METHOD == 'power_fft':
        r = scf.rfft(curr_t*volt_t-curr_a*volt_a)       
        return r
    
    if G.MISFIT_METHOD == 'sko':
        r = (curr_t-curr_a)
        #r2 = r*r
        r2 = np.abs(r)
        return r2
    
    if G.MISFIT_METHOD == 'sko_fft':
        r = (curr_t-curr_a)
        r2 = r*r
        return scf.rfft(r2)
    
    if G.MISFIT_METHOD == 'libivcmp':
        if G.USE_FFT_FILTRATION:
            s = "USE_FFT_FILTRATION=True incompatible with MISFIT_METHOD = 'libivcmp'"
            raise RuntimeError(s)
            
        step_IVCurve = analysis_to_IVCurve()
        res = libivcmp.CompareIvc(target_IVCurve, step_IVCurve, libivcmp.MAX_NUM_POINTS)
        r = np.array([res])
        return r
    
    s = "unknown MISFIT_METHOD = '"+str(G.MISFIT_METHOD)+"'"
    raise RuntimeError(s)

    
    
       
#### ФУНКЦИИ РЕШАТЕЛЯ ########################################################
# вектор, обеспечивающий минимум оптимизируемой функции
Xi_result = [0.,0.,0.,0., 0.,0.,0., 0.,0.,0.,0.] 

# текущий найденный минимум оптимизируемой функции 
misfit_result = 0.

# счетчик числа вызовов функции оптимизатором
FitterCount = 0

# начать работать с решателем - сбросить счетчики, вывести инфо
def init_fitter():
    #print('\ninit_fitter : ')
    global FitterCount
    FitterCount = 0
    

# завершить работу с решателем - вывести инфо
def report_fitter():
    print('\nreport_fitter : ')
    print('FitterCount = '+str(FitterCount))
    print('Xi_result = '+str(Xi_result))
    print('misfit = '+str(misfit_result))
    

BestMisfitCount = 0
# функция вызывается оптимизатором
def fitter_subroutine(Xargs):
    global Xi_result,misfit_result,FitterCount,BestMisfitCount
    FitterCount += 1
    
    xi = Xi_unroll(Xargs)
    
    generate_circuitFile_by_values(xi)
    process_circuitFile()
    
    Mcore = analysis_misfit_core()
    Mscalar = np.sum(np.abs(Mcore))
    Mresult = Mcore
    
    if G.MISFIT_KIND=='minimize': # используем оптимизатор скалярной функции
        Mresult = Mscalar
    
    
    #print("fCount="+str(FitterCount)+', misfit='+str(Mscalar)+', Xargs='+str(Xargs))
    # первый запуск
    if FitterCount<=1:
        Xi_result = xi.copy()
        misfit_result = Mscalar
        BestMisfitCount = 0
        #print("fCount="+str(FitterCount)+', mCount='+str(BestMisfitCount)+', misfit='+str(Mscalar)+', Xargs='+str(Xargs))
            
    # лучший случай
    if Mscalar<misfit_result:
        Xi_result = xi.copy()
        misfit_result = Mscalar
        BestMisfitCount += 1
        #print("fCount="+str(FitterCount)+', mCount='+str(BestMisfitCount)+', misfit='+str(Mscalar)+', Xargs='+str(Xargs))
        
    return Mresult # возвращаем вектор или скаляр


# запустить автоподбор
def run_fitter(result_cir_file_name='',result_csv_file_name=''):       
    if G.MISFIT_KIND == 'minimize': 
        return run_fitter_minimize(result_cir_file_name='',result_csv_file_name='')  
    if G.MISFIT_KIND == 'least_square':        
        return run_fitter_sqleast(result_cir_file_name='',result_csv_file_name='') 
         
    s = "unknown MISFIT_KIND = '"+str(G.MISFIT_KIND)+"'"
    raise RuntimeError(s)
    

# запустить автоподбор - сравнение по набору точек         
def run_fitter_sqleast(result_cir_file_name='',result_csv_file_name=''):    
    global Xi_result,FitterCount             
    print('\nrun_fitter\nXinit = ')
    Xargs = Xi_pack(Xi_long)
    #print('Xargs='+str(Xargs))    
    #print('Xi_mask = '+str(Xi_mask))
    
    FitterCount = 0
     
    x_dif = np.array([1e1,1e1,1e1])
    
    min_bnds = Xargs.copy()
    max_bnds = Xargs.copy()
    
    for i in range(len(min_bnds)):
        min_bnds[i] =np.abs(min_bnds[i])
        max_bnds[i] =np.abs(2.*max_bnds[i])
    
    #min_bnds = (-1e-1,-1e-1,-1e-1)
    #max_bnds = (1e4,1e4,1e4)
    
    bnds = np.array((min_bnds,max_bnds))
    
    #resX = spo.least_squares(fitter_subroutine_sqleast,Xargs,bounds=bnds,xtol=None,ftol=G.TOLERANCE,gtol=None,diff_step=x_dif,max_nfev=1000)
    resX = spo.least_squares(fitter_subroutine,Xargs,jac='3-point',bounds=bnds,method='dogbox',xtol=None,ftol=G.TOLERANCE,gtol=None,diff_step=x_dif,tr_solver='lsmr',max_nfev=1000)
    
    
    print(resX.message)
      
    X = resX.x
    Xi_result = Xi_unroll(X)
    print('FitterCount='+str(FitterCount))
    # вызываем с результатом оптимизации, ибо предыдущий вызов может быть неоптимальным
    generate_circuitFile_by_values(Xi_result)
    process_circuitFile()
              
    if(not result_csv_file_name==''):
        spice.SaveFile(analysis, result_csv_file_name)
    if(not result_cir_file_name==''):          
        generate_circuitFile_by_values(resX.x,result_cir_file_name)
    
    return True


# запустить автоподбор - сравнение по сумме отклонений точек         
def run_fitter_minimize(result_cir_file_name='',result_csv_file_name=''):                 
    # print('\nrun_fitter\nXinit = ')
    Xargs = Xi_pack(Xi_long)
    
    for i in range(0,len(Xargs)):
        Xargs[i] = 0.
       
    resX = spo.minimize(fitter_subroutine,Xargs,method=G.FITTER_METHOD,tol=G.TOLERANCE,options={'maxiter':1000})   
    # вызываем с результатом оптимизации, ибо предыдущий вызов может быть неоптимальным
    # generate_circuitFile_by_values(Xi_result)
    # process_circuitFile()  
              
    if(not result_csv_file_name==''):
        spice.SaveFile(analysis, result_csv_file_name)
    if(not result_cir_file_name==''):          
        generate_circuitFile_by_values(resX.x,result_cir_file_name)
    
    return True
    
### элементарная схема ###
##############################################################################
def Sch_init():
    sch = {}
    sch['R1'] = G.HUGE_R
    sch['C1'] = (G.NONE_C) # [1]
    sch['_R_C1'] = G.NULL_R 
    sch['_R_D1'] = G.HUGE_R
    
    sch['R2'] = G.HUGE_R 
    sch['C2'] = (G.NONE_C) # [5]
    sch['_R_C2'] = G.NULL_R
    
    sch['R3'] = G.HUGE_R 
    sch['C3'] = (G.NONE_C) #[8]
    sch['_R_C3'] = G.NULL_R 
    sch['_R_D3'] = G.HUGE_R         
    return sch


def Sch_set_switch(sch,key,using_sw):
    if using_sw:
        r_sw = G.HUGE_R
    else:
        r_sw = G.NULL_R
        
    if key=='C1': sch['_R_C1']= r_sw
    if key=='C2': sch['_R_C2']= r_sw
    if key=='C3': sch['_R_C3']= r_sw
    if key=='D1': sch['_R_D1']= r_sw
    if key=='D3': sch['_R_D3']= r_sw
    

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


def Sch_init_by_approximation(sch):
    # R1,C1 по максимуму тока и напряжения
    I_max = np.amax(target_VCurrent)
    max_cur_ind = np.argmax(target_VCurrent)
    max_volt_ind = np.argmax(target_input_dummy)
    
    phase_1 = 2.*np.pi*(max_cur_ind-max_volt_ind)/(G.MAX_NUM_POINTS)
    
    V_max = target_input_dummy[max_cur_ind]-G.DIODE_VOLTAGE
    Z1 = V_div_I(V_max,I_max)
    Z1_ = V_div_I(V_max+G.DIODE_VOLTAGE,I_max)
    R1 = Z1*np.cos(phase_1)
    C1 = R_to_C(Z1*np.sin(phase_1))
    
    sch['R1'] = np.abs(R1)
    if np.abs(phase_1*180./np.pi)>1.:
        sch['_R_C1'] = G.HUGE_R
        sch['C1'] = C1
    
    # print('R1='+str(R1)+', C1='+str(C1)+' Z1='+str(Z1)+' ,Z1_='+str(Z1_))  
    # print('phase_1='+str(phase_1*180./np.pi))
    I_min = np.amin(target_VCurrent)
    
    # R3,C3 по минимуму тока и напряжения
    min_cur_ind = np.argmin(target_VCurrent)
    min_volt_ind = np.argmin(target_input_dummy)
    
    V_min = target_input_dummy[min_cur_ind]+G.DIODE_VOLTAGE
    Z3 = V_div_I(V_min,I_min)
    Z3_ = V_div_I(V_min-G.DIODE_VOLTAGE,I_min)
    
    phase_3 = 2.*np.pi*(min_cur_ind-min_volt_ind)/(G.MAX_NUM_POINTS)
    
    R3 = Z3*np.cos(phase_3)
    C3 = R_to_C(Z3*np.sin(phase_3))
    # print('R3='+str(R3)+', C3='+str(C3)+', Z3='+str(Z3)+', Z3_='+str(Z3_))
    # print('phase_3='+str(phase_3*180./np.pi))
    
    sch['R3'] = np.abs(R3)
    if np.abs(phase_3*180./np.pi)>1.:
        sch['_R_C3'] = G.HUGE_R
        sch['C3'] = C3
        
    R2 = R1+R3  # переделать на окрестности малого сигнала
    C2 = C1+C3  # переделать на сдвиг фаз в окрестности ноля
    sch['R2'] = R2
    if C_to_R(C2)>0.01*R2:
        sch['_R_C3'] = G.HUGE_R
        sch['C2'] = C2
        
    
   
    
def init_target_by_Sch(sch):
    generate_circuitFile_by_values(Sch_get_Xi(sch))
    init_target_by_circuitFile()    
    
        
############################################################################# 
def Session_create(start_sch):
    s = {}
    s['start_sch'] = start_sch
    return s


def Session_run(session):
    sch = session['start_sch']
    var_list = session['Xi_variable']
    init_fitter()
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
    r1_mask = 1+8+16 # R1+C1+D1
    r2_mask = 2+128 # R2+C2
    r3_mask = 3+32+64 #R3+C3+D3
    
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
        sch['R1'] = G.HUGE_R
    else:
        var_list +=['R1']
        
    if swcode & 2: # ветка 2 
        #print('dnp branch2') 
        sch['R2'] = G.HUGE_R
    else:
        var_list += ['R2']
    
    if swcode & 4: # ветка 3  
        #print('dnp branch3')
        sch['R3'] = G.HUGE_R
    else:
        var_list += ['R3']
        
    if swcode & 8: # C1  
        #print('dnp C1')
        sch['_R_C1'] = G.NULL_R
    else:
        sch['_R_C1'] = G.HUGE_R
        var_list += ['C1']
        
    if swcode & 16: # D1 
        #print('dnp D1') 
        sch['_R_D1'] = G.NULL_R
    else:
        sch['_R_D1'] = G.HUGE_R
        
    
    if swcode & 32: # C3 
        #print('dnp C3')
        sch['_R_C3'] = G.NULL_R
    else:
        sch['_R_C3'] = G.HUGE_R
        var_list += ['C3']
        
    if swcode & 64: # D3  
        #print('dnp D3')
        sch['_R_D3'] = G.NULL_R
    else:
        sch['_R_D3'] = G.HUGE_R
        
    if swcode & 128: # C2
        #print('dnp C2')
        sch['_R_C2'] = G.NULL_R
    else:
        sch['_R_C2'] = G.HUGE_R
        var_list += ['C2']
    
    session['Xi_variable'] = var_list  
    
    
def Session_processAll():
    ses_result = None
    best_misfit = 1e9
    
    for swcode in range(255):
        if not is_valid_switchers(swcode):
            continue
        
        sch0 = Sch_init()
        Sch_init_by_approximation(sch0)
        ses = Session_create(sch0)
        Session_set_switchers(ses,swcode)
        Session_run(ses)
        if(ses['misfit']<best_misfit):
            best_misfit = ses['misfit']
            ses_result = Session_create(ses['result_sch'])
            Session_set_switchers(ses_result,swcode)
            print('\nswcode = '+str(swcode))
            print('misfit = '+str(best_misfit))
            print(ses['result_sch'])           
            analysis_plot()
                    
    # Session_run(ses_result)
    # analysis_plot()
    
    
#############################################################################


def test1():
    sch = Sch_init()
    sch['R1'] = 1e-2
    sch['C1'] = 1e-5
    sch['_R_C1'] = G.HUGE_R
    sch['R2'] = 1e2
    # sch['_R_C2'] = G.HUGE_R
    # sch['C2'] = 1e-5
    # sch['C3'] = 1.3e-6
    # sch['_R_C3'] = G.HUGE_R
    sch['R3'] = 1e3
    init_target_by_Sch(sch)
    Session_processAll()


def test2():
    sch = Sch_init()
    sch['R1'] = 1e2
    sch['C1'] = 1e-5
    sch['_R_C1'] = G.HUGE_R
    # sch['R2'] = 1e2
    # sch['_R_C2'] = G.HUGE_R
    # sch['C2'] = 1e-5
    # sch['C3'] = 1.3e-6
    # sch['_R_C3'] = G.HUGE_R
    sch['R3'] = 1e3
    init_target_by_Sch(sch)
    Session_processAll()

def test3():
    sch = Sch_init()
    # sch['R1'] = 1e-1
    # sch['C1'] = 1e-5
    # sch['_R_C1'] = G.HUGE_R
    sch['R2'] = 1e2
    sch['_R_C2'] = G.HUGE_R
    sch['C2'] = 1e-7
    # sch['C3'] = 1.3e-6
    # sch['_R_C3'] = G.HUGE_R
    # sch['R3'] = 1e-3
    init_target_by_Sch(sch)
    Session_processAll()
    

def test4():
    sch = Sch_init()
    # sch['R1'] = 1e2
    # sch['C1'] = 1e-5
    # sch['_R_C1'] = G.HUGE_R
    sch['R2'] = 1e2
    sch['_R_C2'] = G.HUGE_R
    sch['C2'] = 1e-7
    # sch['C3'] = 1.3e-6
    # sch['_R_C3'] = G.HUGE_R
    sch['R3'] = 1e3
    init_target_by_Sch(sch)
    Session_processAll()
    
def test_data(csv_data):
    init_target_from_csvFile(csv_data)
    Session_processAll()
      
    
def main():
    set_circuit_template(G.CIR_TEMPLATE)  
    test1()
    test2()
    test3()
    test4()
    test_data('test_data1.csv')
    test_data('test_data2.csv')
    test_data('test_data3.csv')
    
    
    
if __name__=='__main__':
    main()
      
##############################################################################
  
