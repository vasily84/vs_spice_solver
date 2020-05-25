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
# 
import vs_solver_test as my_test

if G.USE_LIBIVCMP:
    G.TOLERANCE = 1e-15
    G.MAX_NUM_POINTS = libivcmp.MAX_NUM_POINTS


# результат последненго моделирования в PySpice
analysis = None

# целевая кривая с током. Та, которую мы подбираем.
target_VCurrent = None
target_input_dummy = None
# заряд, ушедщий в схему.
target_Q = None

# целевая кривая с током для сравнения в библиотеке libivcmp
target_IVCurve = None

# кривая с током для первого приближения
firstStep_VCurrent = None
firstStep_input_dummy = None

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

# Маска оптимизируемых параметров - список булевого типа, например -
# Xi_long = [10.0, x1, 300., x3]
# Xi_mask = [False,True,False,True] -> X_short = [x1,x3]
Xi_mask = [False,False,False,False, False,False,False, False,False,False,False]


#### ФУНКЦИИ ДЛЯ ШАБЛОНА, ЦЕЛЕВОЙ МОДЕЛИ И МАСКИ ПАРАМЕТРОВ ##################
def Xi_unroll(x_short):
    XL = Xi_long
    
    j = 0
    for i in range(0,len(Xi_mask)):
        if Xi_mask[i]:
            XL[i] = x_short[j]
            j += 1
            
    return XL

def Xi_pack():
    xi = []
       
    for i in range(0,len(Xi_mask)):
        if Xi_mask[i]:
            xi += [Xi_long[i]]
                       
    return xi
    


# установить все известные номиналы и маску оптимизации
def set_circuit_nominals_and_mask(nominals,mask):
    global Xi_long,Xi_mask
    # print('set_circuit_nominals_and_mask')
    # print('nominals = '+str(nominals))
    # print('mask = '+str(mask))
    Xi_long = nominals
    # print('Xi_long')
    # print(Xi_long)
    Xi_mask = mask
    
    
# установить файл шаблона схемы. Считать его в набор строк.
def set_circuit_template(fileName):
    global circuit_TempateStrings
    circuit_TempateStrings = []
    templateF = open(fileName)
    
    for tStr in templateF:
        circuit_TempateStrings += [tStr]
        #
    templateF.close()
   
    
# инициализировать целевую модель, промоделироваа файл схемы
def init_target_by_circuitFile(fileName = circuit_SessionFileName):
    global target_VCurrent, target_input_dummy, target_IVCurve,target_Q
    process_circuitFile(fileName)
    target_VCurrent = analysis.VCurrent
    target_input_dummy = analysis.input_dummy
    
    if G.USE_LIBIVCMP:
        target_IVCurve = analysis_to_IVCurve()
        
    target_Q = np.zeros_like(target_VCurrent)
    Q = 0.
    for i in range(0,len(target_VCurrent)):
        Q += target_VCurrent[i]
        target_Q[i] = Q
        
           
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
                #print(target_input_dummy)
            if i==2:
                target_VCurrent = np.array(row,dtype=float) 
                #print(target_VCurrent)
            i += 1
            
    target_Q = np.zeros_like(target_VCurrent)
    Q = 0.
    for i in range(0,len(target_VCurrent)):
        Q += target_VCurrent[i]
        target_Q[i] = Q
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


# в наборе строк шаблона схемы сделать замену {} на значения 
# варьирования Xi_values, сохранить заданным с именем
def generate_circuitFile_by_values( Xi_values, fileCircuit = circuit_SessionFileName):
    newF = open(fileCircuit, 'w')
    i = 0
    for tStr in circuit_TempateStrings:
        cStr = tStr
        if tStr.find('{}')>=0:
            cStr = tStr.format(str(np.abs(Xi_values[i]))) # ставим абсолютные значения номиналов
            i +=1
            
        newF.write(cStr)
    newF.close()

        
# промоделировать файл схемы
def process_circuitFile(fileName=circuit_SessionFileName,csvName=''):
    global analysis
    del(analysis)
    
    circuit = spice.LoadFile(fileName)
    input_data = spice.Init_Data(G.INIT_F, G.INIT_V, G.INIT_Rcs,G.INIT_SNR )
    #analysis = spice.CreateCVC1(circuit, input_data, G.MAX_NUM_POINTS, "input", G.INIT_CYCLE)
    analysis = spice.CreateCVC1(circuit, input_data, G.MAX_NUM_POINTS, "input", G.INIT_CYCLE)
    if(not csvName==''):
        spice.SaveFile(analysis, csvName)
 

# последний анализ перевести в форму, пригодную для сравнения в libivcmp
def analysis_to_IVCurve():
    iv_curve = libivcmp.IvCurve()
    for i in range(G.MAX_NUM_POINTS):
        iv_curve.voltages[i] = c_double(analysis.VCurrent[i])
        iv_curve.currents[i] = c_double(analysis.input_dummy[i])
        
    libivcmp.SetMinVC(0, 0)
    return iv_curve
    
    
# вывести на график результат моделирования
def analysis_plot(title='',pngName=''):
    return
    figure1 = plt.figure(1, (20, 10))
    plt.grid()
    
    # первый шаг итерации подбора ВАХ
    # if not firstStep_input_dummy is None:
    #     plt.plot(firstStep_input_dummy, firstStep_VCurrent,color='yellow')   
    # целевая ВАХ
    #if not target_input_dummy is None:
    plt.plot(target_input_dummy, target_VCurrent,color='red')
    # ВАХ результат подбора
    #if not analysis is None:
    plt.plot(analysis.input_dummy, analysis.VCurrent,color='blue')  
        
    if (not title==''):
       plt.title(title)       
    plt.xlabel('Напряжение [В]')
    plt.ylabel('Сила тока [А]')
    if(not pngName==''):
        plt.savefig(pngName)
        
    plt.show()
    
def analysis_to_Q():
    Q = 0.
    Q_curve = np.zeros_like(analysis.input_dummy)
    for i in range(0,len(analysis.input_dummy)):              
        Q += analysis.VCurrent[i]                
        Q_curve[i] = Q 
        
    return Q_curve
    

# вывести на график результат моделирования
def analysis_plotQ():
    #return 
    figure1 = plt.figure(1, (20, 10))
    plt.grid()
    Q = analysis_to_Q()

    # целевая ВАХ
    plt.plot(target_VCurrent, color='red')
    plt.plot(analysis.VCurrent,color='blue')
    plt.plot(analysis.VCurrent-target_VCurrent,color='black')
              
    plt.show()
                           

#### ФУНКЦИИ СРАВНЕНИЯ ВАХ ###################################################
    
# вычислить несовпадение последнего анализа и целевой функции.
def analysis_misfit_by_Q():
    # разница зарядов
    sko = 0.
    Q = analysis_to_Q()
    for i in range(0,len(target_Q)):
       dQ = target_Q[i]-Q[i]
       sko += np.abs(dQ)
        
    return sko
     
 
def analysis_misfit_by_sko():  
    return analysis_misfit_by_Q()
    # разница токов
    Isub = target_VCurrent-analysis.VCurrent
    ssum = 0.
    # todo - переделать на библиотечный вызов
    for i in range(0,len(analysis.input_dummy)-1):
        # изменение напряжения
        dV = analysis.input_dummy[i+1]-analysis.input_dummy[i] 
        # средний ток
        I = (Isub[i+1]+Isub[i])/2. 
        ssum += np.abs(I*dV)
              
    return ssum


# вычислить несовпадение последнего анализа и целевой функции. 
def analysis_misfit_by_sko_old():  
    # так сравниваем мощности в точках
    s = (target_VCurrent-analysis.VCurrent)*analysis.input_dummy    
    #s = target_VCurrent*target_input_dummy-analysis.VCurrent*analysis.input_dummy    
    #s2 = s*s # работает наилучшим образом
    s2 = np.abs(s) 
    s_sum = np.sum(s2)
    return s_sum


# вычислить несовпадение последнего анализа и целевой функции. 
# использует модуль libivcmp
def analysis_misfit_by_libivcmp():
    step_IVCurve = analysis_to_IVCurve()
    res = libivcmp.CompareIvc(target_IVCurve, step_IVCurve, libivcmp.MAX_NUM_POINTS)
    return res

analysis_misfit = analysis_misfit_by_sko
# установить используемую функция вычисления несовпадения
def restore_analysis_misfit_function():
    global analysis_misfit
    if G.USE_LIBIVCMP:
        analysis_misfit = analysis_misfit_by_libivcmp
    else:
        analysis_misfit = analysis_misfit_by_sko
        
restore_analysis_misfit_function()
    
    
# функция анализа несовпадений ВАХ по фрагментам:

# вычислить несовпадение последнего анализа и целевой функции. 
# Смотрим только положительную полуволну
def analysis_misfit_by_positiveVoltage():
    s = target_VCurrent-analysis.VCurrent  
    s2 = np.abs(s) 
    s_sum = 0.
    for i in range(0,len(target_input_dummy)):
        if target_input_dummy[i]>0.:
            s_sum += s2[i]
          
    return s_sum


# вычислить несовпадение последнего анализа и целевой функции. 
# Смотрим только отрицательную полуволну
def analysis_misfit_by_negativeVoltage():
    s = target_VCurrent-analysis.VCurrent  
    s2 = np.abs(s) 
    s_sum = 0.
    for i in range(0,len(target_input_dummy)):
        if target_input_dummy[i]<0:
            s_sum += s2[i]
          
    return s_sum

       
#### ФУНКЦИИ РЕШАТЕЛЯ ########################################################
# вектор, обеспечивающий минимум оптимизируемой функции
Xi_result = [0.,0.,0.,0., 0.,0.,0., 0.,0.,0.,0.] 
# текущий найденный минимум оптимизируемой функции 
misfit_result = 0.
# счетчик числа вызовов функции оптимизатором
FitterCount = 0

# начать работать с решателем - сбросить счетчики, вывести инфо
def init_fitter():
    global FitterCount,analysis_misfit_best
    FitterCount = 0
    #print('\ninit_fitter : ')
    

# завершить работу с решателем - вывести инфо
def report_fitter():
    # print('\nreport_fitter : ')
    # print('FitterCount = '+str(FitterCount))
    # print('Xi_result = '+str(Xi_result))
    # print('misfit = '+str(misfit_result))
    pass
    

# функция вызывается оптимизатором
def fitter_subroutine_sqleast(Xargs):
    global FitterCount,misfit_result,Xi_result,firstStep_VCurrent,firstStep_input_dummy  
    Xargs =1000.*Xargs # всё в килоомах
    xi = Xi_unroll(Xargs)
    #print('xi='+str(xi))
    generate_circuitFile_by_values(xi)
    process_circuitFile()
    
    res = target_Q
    Q = analysis_to_Q()
    A = target_VCurrent#*target_input_dummy
    B = analysis.VCurrent#*analysis.input_dummy
    
    # for i in range(1,len(target_input_dummy)):
    #     A = target_VCurrent[i]*(target_input_dummy[i]-target_input_dummy[i-1])
    #     B = analysis.VCurrent[i]*(analysis.input_dummy[i]-analysis.input_dummy[i-1])
    
    
    FitterCount += 1
    R = (B-A)/G.MAX_NUM_POINTS
    #R = Q-res
    return (R) # возвращаем вектор разностей значений заряда


# функция вызывается оптимизатором
def fitter_subroutine_minimize(Xargs):
    global FitterCount,misfit_result,Xi_result,firstStep_VCurrent,firstStep_input_dummy  
    xi = Xi_unroll(Xargs)
    generate_circuitFile_by_values(xi)
    process_circuitFile()
    
    M = analysis_misfit()
    if G.USE_LIBIVCMP:
        Minv = analysis_misfit_by_libivcmp()
    else:
        Minv = analysis_misfit_by_sko() 
        
    FitterCount += 1
    #print('FitterCount = '+str(FitterCount))
    
    # первый вызов этой функции
    if FitterCount==1:  
        Xi_result = xi
        misfit_result = Minv       
        firstStep_VCurrent = np.copy(analysis.VCurrent)
        firstStep_input_dummy = np.copy(analysis.input_dummy)
        return M
    
    # нашли минимум, лучший чем предыдущие, запоминаем его
    if Minv<misfit_result:      
        Xi_result = xi
        misfit_result = Minv
        
    return M

def run_fitter(result_cir_file_name='',result_csv_file_name=''):       
    return run_fitter_sqleast(result_cir_file_name='',result_csv_file_name='')          

# запустить автоподбор         
def run_fitter_sqleast(result_cir_file_name='',result_csv_file_name=''):    
    global Xi_result,FitterCount             
    print('\nrun_fitter\nXinit = ')
    Xargs = Xi_pack()
    # print('Xargs='+str(Xargs))
    # print(Xargs)
    # print('Xi_mask = ')
    # print(Xi_mask)
    FitterCount = 0
    Xargs[0] = 1e-3*Xargs[0]
    Xargs[1] = 1e-3*Xargs[1]
    Xargs[2] = 1e-3*Xargs[2]
    
    #x_scale = np.array([1e-3,1e-3,1e-3])
    x_dif = np.array([5e-2,5e-2,5e-2])
    bnds = np.array(([-1e2,-1e2,-1e2],[1e8,1e8,1e8]))
    #resX = spo.least_squares(fitter_subroutine_sqleast,Xargs,method='lm',xtol=1e-15)
    resX = spo.least_squares(fitter_subroutine_sqleast,Xargs,bounds=bnds,xtol=None,ftol=None,gtol=1e-5,loss='soft_l1',diff_step=x_dif,max_nfev=1000)
    #print(resX)
    print(resX.message)
    x = resX.x
    #x = x*0.5
    #resX = spo.least_squares(fitter_subroutine_sqleast,x,bounds=bnds,xtol=None,ftol=None,gtol=1e-5,loss='soft_l1',diff_step=x_dif,max_nfev=1000)
    
    X =1000.*resX.x
    Xi_result = Xi_unroll(X)
    #print('resX.x='+str(resX.x))
    print('FitterCount='+str(FitterCount))
    # вызываем с результатом оптимизации, ибо предыдущий вызов может быть неоптимальным
    generate_circuitFile_by_values(Xi_result)
    process_circuitFile()
    #report_fitter()
              
    if(not result_csv_file_name==''):
        spice.SaveFile(analysis, result_csv_file_name)
    if(not result_cir_file_name==''):          
        generate_circuitFile_by_values(resX.x,result_cir_file_name)
    
    return True

# запустить автоподбор         
def run_fitter_minimize(result_cir_file_name='',result_csv_file_name=''):                 
    print('\nrun_fitter\nXinit = ')
    Xargs = Xi_pack()
    print(Xargs)
    print('Xi_mask = ')
    print(Xi_mask)
    
    resX = spo.minimize(fitter_subroutine_minimize,Xargs,method=G.FITTER_METHOD,tol=G.TOLERANCE,options={'maxiter':100})
    #print(resX.message)
    
    # вызываем с результатом оптимизации, ибо предыдущий вызов может быть неоптимальным
    generate_circuitFile_by_values(Xi_result)
    process_circuitFile()
    report_fitter()
              
    if(not result_csv_file_name==''):
        spice.SaveFile(analysis, result_csv_file_name)
    if(not result_cir_file_name==''):          
        generate_circuitFile_by_values(resX.x,result_cir_file_name)
    
    return True
    


##############################################################################

#### класс общей схемы обратной задачи ######################################
class CGeneralCircuit():
    def __init__(self):
        # все ветки отключены
        self.R1 = G.HUGE_R # [0]
        self.C1 = G.NONE_C # [1]
        self.R_C1 = G.NULL_R # [2]
        self.R_D1 = G.HUGE_R # [3]
        
        self.R2 = G.HUGE_R # [4]
        self.C2 = G.NONE_C # [5]
        self.R_C2 = G.NULL_R # [6]
        
        self.R3 = G.HUGE_R # [7]
        self.C3 = G.NONE_C # [8]
        self.R_C3 = G.NULL_R # [9]
        self.R_D3 = G.HUGE_R # [10]
        
        # считаем, что все ветки присутствуют
        self.node_1 = True
        self.node_2 = True
        self.node_3 = True
        
        self.Text = ''
        self.vmask = [False,False,False,False, False,False,False, False,False,False,False]
        
    def reset_values(self):
        # все ветки отключены
        self.R1 = G.HUGE_R # [0]
        self.C1 = G.NONE_C # [1]
        self.R_C1 = G.NULL_R # [2]
        self.R_D1 = G.HUGE_R # [3]
        
        self.R2 = G.HUGE_R # [4]
        self.C2 = G.NONE_C # [5]
        self.R_C2 = G.NULL_R # [6]
        
        self.R3 = G.HUGE_R # [7]
        self.C3 = G.NONE_C # [8]
        self.R_C3 = G.NULL_R # [9]
        self.R_D3 = G.HUGE_R # [10]
        
        self.node_1 = True
        self.node_2 = True
        self.node_3 = True
        
        
    def __repr__(self):
        s1 = ' R1='+str(self.R1)+' C1='+str(self.C1)+' R_C1='+str(self.R_C1)+' R_D1='+str(self.R_D1)
        s2 = ' R2='+str(self.R2)+' C2='+str(self.C2)+' R_C2='+str(self.R_C1)
        s3 = ' R3='+str(self.R3)+' C3='+str(self.C3)+' R_C3='+str(self.R_C3)+' R_D3='+str(self.R_D3)
        return (s1+s2+s3)
    
    # установить номиналы в вектор для оптимизации 
    def setup(self):
        set_circuit_template('general.cir_t')                   
        set_circuit_nominals_and_mask(self.get_Xi(),self.vmask)
       
    def load_result(self):
        # print('\nload_result:')
        # print(Xi_result)
        if self.vmask[0]: self.R1 = Xi_result[0]
        if self.vmask[1]: self.C1 = Xi_result[1]
        if self.vmask[2]: self.R_C1 = Xi_result[2]
        if self.vmask[3]: self.R_D1 = Xi_result[3]
        
        if self.vmask[4]: self.R2 = Xi_result[4]
        if self.vmask[5]: self.C2 = Xi_result[5]
        if self.vmask[6]: self.R_C2 = Xi_result[6]
        
        if self.vmask[7]: self.R3 = Xi_result[7]
        if self.vmask[8]: self.C3 = Xi_result[8]
        if self.vmask[9]: self.R_C3 = Xi_result[9]
        if self.vmask[10]: self.R_D3 = Xi_result[10]
        
       
    def reset_vmask(self):
        self.vmask = [False,False,False,False, False,False,False, False,False,False,False]
        
        
    # сформировать вектор из значений.
    def get_Xi(self):
        xi = []
        xi += [self.R1]
        xi += [self.C1]
        xi += [self.R_C1]
        xi += [self.R_D1]
        
        xi += [self.R2]
        xi += [self.C2]
        xi += [self.R_C2]
        
        xi += [self.R3]
        xi += [self.C3]
        xi += [self.R_C3]
        xi += [self.R_D3]
        
        return xi
    
    def init_target_random(self):
        self.reset_values()
        self.R1 = 1e-5      
        #self.R2 = 1e3*(random.random())            
        self.R3 = 1.
            
        self.setup()        
        f = 'target1.cir'
        generate_circuitFile_by_values(self.get_Xi(), f)
        init_target_by_circuitFile(f)
        
        
    def init_target1(self):
        self.reset_values()
        #self.R1 = 100 
        #self.C1 = 1e-6
        #self.R_C1 = G.NULL_R
        
        #self.R3 = 100 
        #self.C3 = 1e-6
        #self.R_C3 = G.NULL_R
                   
        self.R2 = 3000
        self.setup()
        
        generate_circuitFile_by_values(self.get_Xi())
        init_target_by_circuitFile()
    
    def init_target2(self):
        self.reset_values()    
        self.R1 = 2.
        self.R3 = 1000   
        #self.R_C3 = G.HUGE_R       
        #self.C3= 1e-7     
        self.R2 = 400                             
        self.setup()
        
        generate_circuitFile_by_values(self.get_Xi())
        init_target_by_circuitFile()
    
    def init_target3(self):
        self.reset_values()  
        self.R1 = 1.       
        self.R2 = 100
        #self.R_C2 = G.HUGE_R
        #self.C2 = 1e-6
        self.R3 = 1200
        #self.C2 = 1e-6  
        #self.R_C2 = G.HUGE_R                                     
        self.setup()
        
        generate_circuitFile_by_values(self.get_Xi())
        init_target_by_circuitFile()
        
        
    # определить сопротивление на положительной полуволне
    def get_target_positiveVoltage_Resistance(self):       
        Rsumm =0.
        Rcount = 0
        
        for i in range(0,len(target_input_dummy)):
            if target_input_dummy[i]>G.SMALL_VOLTAGE:
                I= target_VCurrent[i]
                V= target_input_dummy[i]      
                try:
                    R = np.abs(V/I)
                except ArithmeticError:
                    R = G.HUGE_R
                Rsumm += R
                Rcount += 1
                
        try:
            R = Rsumm/Rcount
        except ArithmeticError:
            R = G.HUGE_R
            
        return R
    
    # определить сопротивление на отрицательной полуволне
    def get_target_negativeVoltage_Resistance(self):
        Rsumm =0.
        Rcount = 0
        
        for i in range(0,len(target_input_dummy)):
            if target_input_dummy[i]<-1.*G.SMALL_VOLTAGE:
                I= target_VCurrent[i]
                V= target_input_dummy[i]      
                try:
                    R = np.abs(V/I)
                except ArithmeticError:
                    R = G.HUGE_R
                Rsumm += R
                Rcount += 1
                
        try:
            R = Rsumm/Rcount
        except ArithmeticError:
            R = G.HUGE_R
            
        return R


    # определить сопротивление на положительной полуволне, как если бы она была
    # включена через диод. Возвращает положительные и отрицательные значения
    def get_target_positiveVoltage_Resistance_via_Diode(self):
        Rsumm =0.
        Rcount = 0
        
        for i in range(0,len(target_input_dummy)):
            if target_input_dummy[i]>G.DIODE_VOLTAGE:
                I= target_VCurrent[i]
                V= target_input_dummy[i]-G.DIODE_VOLTAGE      
                try:
                    R = np.abs(V/I)
                except ArithmeticError:
                    R = G.HUGE_R
                Rsumm += R
                Rcount += 1
                
        try:
            R = Rsumm/Rcount
        except ArithmeticError:
            R = G.HUGE_R
            
        return R
    
    # определить сопротивление на отрицательной полуволне, как если бы она была
    # включена через диод. Возвращает положительные и отрицательные значения
    def get_target_negativeVoltage_Resistance_via_Diode(self):
        Rsumm =0.
        Rcount = 0
        
        for i in range(0,len(target_input_dummy)):
            if target_input_dummy[i]<-1.*G.DIODE_VOLTAGE:
                I= target_VCurrent[i]
                V= target_input_dummy[i]+G.DIODE_VOLTAGE      
                try:
                    R = np.abs(V/I)
                except ArithmeticError:
                    R = G.HUGE_R
                Rsumm += R
                Rcount += 1
                
        try:
            R = Rsumm/Rcount
        except ArithmeticError:
            R = G.HUGE_R
            
        return R
    
    # определить сопротивление в окрестности нуля
    def get_target_smallVoltage_Resistance(self):
        Rsumm =0.
        Rcount = 0
        
        for i in range(0,len(target_input_dummy)):
            v = np.abs(target_input_dummy[i])
            if v>=G.SMALL_VOLTAGE and v<=G.DIODE_VOLTAGE:
                I= target_VCurrent[i]
                V= target_input_dummy[i]      
                try:
                    R = np.abs(V/I)
                except ArithmeticError:
                    R = G.HUGE_R
                Rsumm += R
                Rcount += 1
                
        try:
            R = Rsumm/Rcount
        except ArithmeticError:
            R = G.HUGE_R
            
        return R
    
    # определить сопротивление по всему диапазону напряжений
    def get_target_total_Resistance(self):
        Rsumm =0.
        Rcount = 0
        
        for i in range(0,len(target_input_dummy)):
            if np.abs(target_input_dummy[i])>G.SMALL_VOLTAGE:
                I= target_VCurrent[i]
                V= target_input_dummy[i]      
                try:
                    R = np.abs(V/I)
                except ArithmeticError:
                    R = G.HUGE_R
                Rsumm += R
                Rcount += 1
                
        try:
            R = Rsumm/Rcount
        except ArithmeticError:
            R = G.HUGE_R
            
        return R
    
    
    # найти стартовые значения для оптимизации и определить присутствие цепи 1,2,3
    def find_initial_resistors_and_detect_nodes(self):
        self.reset_vmask()
        r_plus = self.get_target_positiveVoltage_Resistance()
        r_plus_via_diode = self.get_target_positiveVoltage_Resistance_via_Diode()
        r_minus = self.get_target_negativeVoltage_Resistance()
        r_minus_via_diode = self.get_target_negativeVoltage_Resistance_via_Diode()
        r_small = self.get_target_smallVoltage_Resistance()
        print('\nmeasured values: \n')
        print('r+ = '+str(r_plus))
        print('r- = '+str(r_minus))
        print('r+_diode = '+str(r_plus_via_diode))
        print('r-_diode = '+str(r_minus_via_diode))        
        print('r_small = '+str(r_small))
                                   
        #sigma_2 = 1./r_small # старый вариант         
        sigma_12 = 1./(r_plus_via_diode)
        sigma_23 = 1./(r_minus_via_diode)
        sigma_total = 1./self.get_target_total_Resistance()  
        sigma_2 = sigma_total-sigma_12-sigma_23
                
        print('initial values:\n')
        # ветки 1 и 2 отсутствуют
        if (r_plus>=G.LARGE_R):
            self.node_1 = False
            self.node_2 = False
            self.R3 = r_minus_via_diode
            self.vmask[7] = True
            self.setup()
            print('only node_3, r3 = '+str(self.R3))
            return
            
        # ветки 2 и 3 отсутствуют
        if (r_minus>=G.LARGE_R):
            self.node_3 = False
            self.node_2 = False
            self.R1 = r_plus_via_diode
            self.vmask[0] = True
            print('only node_1, r1 = '+str(self.R1))
            self.setup()
            return
        
        # присутствуют все ветки, тяжелый вычислительный случай.
        self.R3 = 1.#1./(sigma_12-sigma_2)
        self.R2 = 1.#r_small
        self.R1 = 1.#1./(sigma_23-sigma_2)
        self.vmask[0] = True
        self.vmask[4] = True
        self.vmask[7] = True
        
        print('r1='+str(self.R1)+', r2='+str(self.R2)+', r3='+str(self.R3))
        self.setup()
        
            
        
    # запустить подбор значений R1,R2, используется только 
    # положительная полуволна напряжения
    def run_fitter_r1r2(self):
        global analysis_misfit 
                   
        self.reset_vmask()
        self.vmask[0] = True # варьируем R1
        self.vmask[4] = True # варьируем R2
        
        analysis_misfit = analysis_misfit_by_positiveVoltage
        self.setup()
        
        run_fitter()                        
        self.load_result()    
                    
        
        
    # запустить подбор значений R2,R3. используется только отрицательная 
    # полуволна напряжения 
    def run_fitter_r2r3(self):
        global analysis_misfit 
                    
        self.reset_vmask()
        self.vmask[4] = True # R2
        self.vmask[7] = True # R3
        
        analysis_misfit = analysis_misfit_by_negativeVoltage
        self.setup()
        
        run_fitter()
        self.load_result()     
               
        
        
    # запустить подбор значений R1,R2,R3
    def run_fitter_r1r2r3(self):
        global analysis_misfit 
            
        self.reset_vmask()
        self.vmask[0] = True # R1
        self.vmask[4] = True # R2
        self.vmask[7] = True # R3
        
        restore_analysis_misfit_function()
        
        self.setup()
        
        run_fitter()
        self.load_result()    
            
        
    # запустить подбор значений R1
    def run_fitter_r1(self):
        global analysis_misfit 
            
        self.reset_vmask()
        self.vmask[0] = True # R1
          
        analysis_misfit = analysis_misfit_by_positiveVoltage       
        self.setup()
        
        run_fitter()  
        self.load_result()    
        
        
        
    # запустить подбор значений R2
    def run_fitter_r2(self):
        global analysis_misfit 
            
        self.reset_vmask()
        self.vmask[4] = True # R2       
        restore_analysis_misfit_function()   
        self.setup()      
        run_fitter()
        self.load_result()                 
        
    
        
    # запустить подбор значений R3
    def run_fitter_r3(self):
        global analysis_misfit 
            
        self.reset_vmask()        
        self.vmask[7] = True # R3
        
        analysis_misfit = analysis_misfit_by_negativeVoltage
        self.setup()
        
        run_fitter()
        self.load_result()    
       
          

    # запустить подбор значений R1,R3
    def run_fitter_r1r3(self):
        global analysis_misfit 
            
        self.reset_vmask()
        self.vmask[0] = True # R1
        self.vmask[7] = True # R3
        
        analysis_misfit = analysis_misfit_by_sko
        self.setup()
        
        run_fitter()
        self.load_result()   
    
    # запустить подбор значений для с2
    def run_fitter_c2(self):
        global analysis_misfit 
            
        self.reset_vmask()
        self.vmask[5] = True # C2
        
        analysis_misfit = analysis_misfit_by_sko
        self.setup()
        
        run_fitter()
        self.load_result()   
        
        
    # сброс внутренних значений, обнуление для старта    
    def reset_fitter(self):
        self.reset_values()
        self.vmask=[True,True,True,True, True,True,True, True,True,True,True]
        self.setup()
        
    def analysis(self):
        f = 'target3.cir'
        generate_circuitFile_by_values(self.get_Xi(),f)
        
    def run_fitter(self):
        self.reset_fitter()
        self.find_initial_resistors_and_detect_nodes()
        init_fitter()
        #G.INIT_SNR = 100.
        # рассматриваем частные случаи                  
        # только ветка 1
        # if self.node_1==True and self.node_2==False and self.node_3==False:
        #     print('node_1 calculation:')
        #     self.run_fitter_r1()
        #     report_fitter()
        #     analysis_plot(title='РЕЗУЛЬТАТ')            
        #     return
        
        
        # # только ветка 3
        # if self.node_3==True and self.node_1==False and self.node_2==False:
        #     print('node_2 calculation:')
        #     self.run_fitter_r3()
        #     report_fitter()
        #     analysis_plot(title='РЕЗУЛЬТАТ')             
        #     return
        
        # # только ветка 2, нам повезло
        # if self.node_2==True and self.node_1==False and self.node_3==False:
        #     print('node_2 calculation:')
        #     self.run_fitter_r2()
        #     report_fitter()
        #     analysis_plot(title='РЕЗУЛЬТАТ')              
        #     return
        
        
        # все остальные случаи  
                      
        #self.run_fitter_r1r2r3()  
        self.run_fitter_r1r2r3()  
        report_fitter()   
        analysis_plot(title='РЕЗУЛЬТАТ')
        analysis_plotQ()
        
            
        
    def run_me(self):           
        # self.init_target_random() # тест ветки 1
        # self.reset_values()    
        # self.run_fitter()
        #G.restore_INIT_param()   
        self.init_target1() # тест ветки 1
        self.reset_values()    
        self.run_fitter()
        G.restore_INIT_param()
         
        self.init_target2()       
        self.reset_values()
        self.run_fitter()
        G.restore_INIT_param()
       
        # self.init_target3()    
        # analysis_plotQ()
        # analysis_plot()
        # #analysis_plot()
        # self.reset_values()
        # #self.run_fitter()
        
        
#############################################################################    
Sch = CGeneralCircuit()
def main():
    #my_test.test_all()
    #case1 = test_case_1()
    #print(case1)
    
    Sch.run_me()
    
    
if __name__=='__main__':
    main()
      
##############################################################################
  
