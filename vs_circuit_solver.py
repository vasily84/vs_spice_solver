# vs_circuit_solver.py
# номер версии не присвоен
# язык Python
#
# программа подбора значений R,C для вариантов электронной схемы
# исходя из моделирования подобной схемы в ngspice
# поставляется без всякой оптимизации, ибо имеет целью установление методики 
# расчета таких вещей и определения границ применимости этой методики
#
# автор В.Симонов, 21-мая-2020
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
    global Xi_long
    XL = Xi_long
    
    j = 0
    for i in range(0,len(Xi_mask)):
        if Xi_mask[i]:
            XL[i] = x_short[j]
            j += 1
            
    return XL


# установить все известные номиналы и маску оптимизации
def set_circuit_nominals_and_mask(nominals,mask):
    global Xi_long,Xi_mask
    Xi_long = nominals
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
    global target_VCurrent, target_input_dummy, target_IVCurve
    process_circuitFile(fileName)
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
                print(target_input_dummy)
            if i==2:
                target_VCurrent = np.array(row,dtype=float) 
                print(target_VCurrent)
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
    circuit = spice.LoadFile(fileName)
    input_data = spice.Init_Data(G.INIT_F, G.INIT_V, G.INIT_Rcs,G.INIT_SNR )
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
    figure1 = plt.figure(1, (20, 10))
    plt.grid()
    plt.plot(target_input_dummy, target_VCurrent,color='red')
    plt.plot(analysis.input_dummy, analysis.VCurrent,color='blue')  
    #plt.plot(firstStep_input_dummy, firstStep_VCurrent,color='yellow')
        
    if (not title==''):
       plt.title(title)       
    plt.xlabel('Напряжение [В]')
    plt.ylabel('Сила тока [А]')
    if(not pngName==''):
        plt.savefig(pngName)
        
    plt.show()
                           

#### ФУНКЦИИ СРАВНЕНИЯ ВАХ ###################################################
    
# вычислить несовпадение последнего анализа и целевой функции. 
# возвращает неотрицательное число    
def analysis_misfit_by_sko():   
    s = target_VCurrent-analysis.VCurrent    
    #s2 = s*s # работает наилучшим образом
    s2 = np.abs(s) 
    s_sum = np.sum(s2)
    return s_sum

# специальная метрика - чем дальше от нуля, тем больший вес имеет
def analysis_misfit_by_sko_v2():   
    s = target_VCurrent-analysis.VCurrent    
    s2 = s*s 
    s_sum = 0.
    for i in range(0,len(analysis.input_dummy)):
        if np.abs(analysis.VCurrent[i])>=G.DIODE_VOLTAGE:
            s_sum += s2[i]
    
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

# счетчик числа вызова функции оптимизатором
ffCount = 0
# функция вызывается оптимизатором
def fitter_subroutine(Xargs):
    global ffCount  
    generate_circuitFile_by_values(Xi_unroll(Xargs))
    process_circuitFile()
    ffCount += 1
    m = analysis_misfit()
    #print('fitter_subroutine Count = '+str(ffCount))
    
    return m


# запустить автоподбор         
def run_fitter(result_cir_file_name='',result_csv_file_name='',Xargs = Xi_long):
    global ffCount, firstStep_VCurrent, firstStep_input_dummy
    ffCount = 0
      
    
    fitter_subroutine(Xargs)
    firstStep_VCurrent = analysis.VCurrent
    firstStep_input_dummy = analysis.input_dummy
    
    resX = spo.minimize(fitter_subroutine,Xargs,method=G.FITTER_METHOD,tol=G.TOLERANCE,options={'maxiter':100})
    
    # вызываем с результатом оптимизации, ибо предыдущий вызов может быть неоптимальным
    fitter_subroutine(resX.x) 
    
    print(resX.message)
    #print('Rtarg = '+str(get_target_positiveVoltage_Resistance()))
    #print('Rtarg_via_diode = '+str(get_target_positiveVoltage_Resistance_via_Diode()))
    
    if resX.success:
        if(not result_csv_file_name==''):
            spice.SaveFile(analysis, result_csv_file_name)
        if(not result_cir_file_name==''):          
            generate_circuitFile_by_values(resX.x,result_cir_file_name,)
    
    return resX
    


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
    
    # установить для моделирования 
    def setup(self):
        set_circuit_template('general.cir_t')       
        set_circuit_nominals_and_mask(self.get_Xi(),self.vmask)
       
        
        
    def reset_vmask(self):
        self.vmask = [False,False,False,False, False,False,False, False,False,False,False]
        
    # # установить случайны номиналы
    # def randomize_all_resistors(self):
    #     self.R1 = 1e6*random.random()
    #     self.R2 = 1e6*random.random()
    #     self.R3 = 1e6*random.random()
        
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
        self.R1 = 1e3*random.random()        
        self.R2 = 1e3*(1+random.random())            
        self.R3 = 2e3*random.random()
            
        self.setup()        
        f = 'target1.cir'
        generate_circuitFile_by_values(self.get_Xi(), f)
        init_target_by_circuitFile(f)
        
        
    def init_target1(self):
        self.reset_values()
        self.R1 = 1e3                                   
        self.setup()
        
        f = 'target1.cir'
        generate_circuitFile_by_values(self.get_Xi(), f)
        init_target_by_circuitFile(f)
    
    def init_target2(self):
        self.reset_values()        
        self.R3 = 1e4                                       
        self.setup()
        
        f = 'target2.cir'
        generate_circuitFile_by_values(self.get_Xi(), f)
        init_target_by_circuitFile(f)
    
    def init_target3(self):
        self.reset_values()  
        self.R1 = 100                
        self.R2 = 100                                             
        self.setup()
        
        f = 'target3.cir'
        generate_circuitFile_by_values(self.get_Xi(), f)
        init_target_by_circuitFile(f)
        
        
    # определить сопротивление на положительной полуволне
    def get_target_positiveVoltage_Resistance(self):
        Isumm = 0.
        Vsumm = 0.
        
        for i in range(0,len(target_input_dummy)):
            if target_input_dummy[i]>0.:
                Isumm += target_VCurrent[i]
                Vsumm += target_input_dummy[i]      
        try:
            R = Vsumm/Isumm
        except ArithmeticError:
            R = G.HUGE_R
            
        return R 
    
    # определить сопротивление на отрицательной полуволне
    def get_target_negativeVoltage_Resistance(self):
        Isumm = 0.
        Vsumm = 0.
        
        for i in range(0,len(target_input_dummy)):
            if target_input_dummy[i]<0.:
                Isumm += target_VCurrent[i]
                Vsumm += target_input_dummy[i]      
        try:
            R = np.abs(Vsumm/Isumm)
        except ArithmeticError:
            R = G.HUGE_R
            
        return R 


    # определить сопротивление на положительной полуволне, как если бы она была
    # включена через диод. Возвращает положительные и отрицательные значения
    def get_target_positiveVoltage_Resistance_via_Diode(self):
        Isumm = 0.
        Vsumm = 0.
       
        for i in range(0,len(target_input_dummy)):
            if target_input_dummy[i]>G.DIODE_VOLTAGE:
                Isumm += target_VCurrent[i]
                Vsumm += (target_input_dummy[i]-G.DIODE_VOLTAGE)          
        try:
            R = Vsumm/Isumm
        except ArithmeticError:
            R = G.HUGE_R
            
        return R 
    
    # определить сопротивление на отрицательной полуволне, как если бы она была
    # включена через диод. Возвращает положительные и отрицательные значения
    def get_target_negativeVoltage_Resistance_via_Diode(self):
        Isumm = 0.
        Vsumm = 0.
       
        for i in range(0,len(target_input_dummy)):
            if target_input_dummy[i]<-1.*G.DIODE_VOLTAGE:
                Isumm += target_VCurrent[i]
                Vsumm += (target_input_dummy[i]+G.DIODE_VOLTAGE) # знак +, ибо само напряжение отрицательное         
        try:
            R = np.abs(Vsumm/Isumm)
        except ArithmeticError:
            R = G.HUGE_R
            
        return R 
    
    # 
    def get_target_smallVoltage_Resistance(self):
        Isumm = 0.
        Vsumm = 0.
       
        for i in range(0,len(target_input_dummy)):
            if np.abs(target_input_dummy[i])<G.SMALL_VOLTAGE:
                Isumm += np.abs(target_VCurrent[i])
                Vsumm += np.abs(target_input_dummy[i])         
        try:
            R = Vsumm/Isumm
        except ArithmeticError:
            R = G.HUGE_R
            
        return R 
    
    #
    def get_target_total_Resistance(self):
        Isumm = np.sum(np.abs(target_VCurrent))
        Vsumm = np.sum(np.abs(target_input_dummy))         
        try:
            R = Vsumm/Isumm
        except ArithmeticError:
            R = G.HUGE_R
            
        return R
    
    
    # найти стартовые значения для оптимизации и определить присутствие цепи 1,2,3
    def find_initial_value_and_detect_nodes(self):
        r_plus = self.get_target_positiveVoltage_Resistance()
        r_plus_via_diode = self.get_target_positiveVoltage_Resistance_via_Diode()
        r_minus = self.get_target_negativeVoltage_Resistance()
        r_minus_via_diode = self.get_target_negativeVoltage_Resistance_via_Diode()
        print('\nmeasured values: \n')
        print('r+ = '+str(r_plus))
        print('r- = '+str(r_minus))
        print('r+_diode = '+str(r_plus_via_diode))
        print('r-_diode = '+str(r_minus_via_diode))
        r_small = self.get_target_smallVoltage_Resistance()
        print('r_small = '+str(r_small))
                                   
        #sigma_2 = 1./r_small
        sigma_total = 1./self.get_target_total_Resistance()        
        sigma_12 = 1./(r_plus_via_diode)
        sigma_23 = 1./(r_minus_via_diode)
        sigma_2 = sigma_total-sigma_12-sigma_23
        r_2 = 1./sigma_2
        
        print('initial values:\n')
        # ветки 1 и 2 отсутствуют
        if (r_plus>=G.LARGE_R)or(r_plus>10.*r_minus):
            self.node_1 = False
            self.node_2 = False
            self.R3 = r_minus_via_diode
            self.setup()
            print('only node_3, r3 = '+str(self.R3))
            return
            
        # ветки 2 и 3 отсутствуют
        if (r_minus>=G.LARGE_R)or(r_minus>10.*r_plus):
            self.node_3 = False
            self.node_2 = False
            self.R1 = r_plus_via_diode
            print('only node_1, r1 = '+str(self.R1))
            self.setup()
            return
        
         # присутствуют все ветки, тяжелый случай.
        self.R3 = 1./(sigma_12-sigma_2)
        self.R2 = r_small
        self.R1 = 1./(sigma_23-sigma_2)
        print('r1='+str(self.R1)+', r2='+str(self.R2)+', r3='+str(self.R3))
        self.setup()
        
            
        
    # запустить подбор значений R1,R2, используется только 
    # положительная полуволна напряжения
    def run_fitter_node12(self):
        global analysis_misfit 
                   
        self.reset_vmask()
        self.vmask[0] = True # варьируем R1
        self.vmask[4] = True # варьируем R2
        
        analysis_misfit = analysis_misfit_by_positiveVoltage
        self.setup()
        
        res = run_fitter()
        
        r1 = res.x[0]
        r2 = res.x[1]
        print('r1 = '+str(r1))
        print('r2 = '+str(r2))
        
        if res.success: 
            self.R1 = r1
            self.R2 = r2
            self.setup()    # сохраняем вектор с успешным подбором
            
        print(self)
        analysis_plot(title='подбор R1,R2\n по положительной полуволне')
        
        
    # запустить подбор значений R2,R3. используется только отрицательная 
    # полуволна напряжения 
    def run_fitter_node23(self):
        global analysis_misfit 
                    
        self.reset_vmask()
        self.vmask[4] = True # R2
        self.vmask[7] = True # R3
        
        analysis_misfit = analysis_misfit_by_negativeVoltage
        self.setup()
        
        res = run_fitter()
        
        r2 = np.abs(res.x[0])
        r3 = np.abs(res.x[1])
        print('r2 = '+str(r2))
        print('r3 = '+str(r3))
        
        if res.success: 
            self.R2 = r2
            self.R3 = r3
            self.setup()    # сохраняем вектор с успешным подбором
            

        print(self)
        analysis_plot(title='подбор R2,R3\n по отрицательной полуволне')
        
    # запустить подбор значений R1,R2,R3
    def run_fitter_node123(self):
        global analysis_misfit 
            
        self.reset_vmask()
        self.vmask[0] = True # R1
        self.vmask[4] = True # R2
        self.vmask[7] = True # R3
        
        restore_analysis_misfit_function()
        
        self.setup()
        
        res = run_fitter()
        
        r1 = np.abs(res.x[0])
        r2 = np.abs(res.x[1])
        r3 = np.abs(res.x[2])
        print('r1 = '+str(r1))
        print('r2 = '+str(r2))
        print('r3 = '+str(r3))
        
        if res.success: 
            self.R1 = r1
            self.R2 = r2
            self.R3 = r3
            self.setup()    # сохраняем вектор с успешным подбором
            

        print(self)
        analysis_plot(title='подбор R1,R2,R3')
        
    # запустить подбор значений R1
    def run_fitter_node1(self):
        global analysis_misfit 
            
        self.reset_vmask()
        self.vmask[0] = True # R1
          
        analysis_misfit = analysis_misfit_by_positiveVoltage       
        self.setup()
        
        res = run_fitter()       
        r1 = np.abs(res.x[0])   
        print('r1 = '+str(r1))
       
        
        if res.success: 
            self.R1 = r1         
            self.setup()    # сохраняем вектор с успешным подбором
            

        print(self)
        analysis_plot(title='подбор R1')
        
        
    # запустить подбор значений R2
    def run_fitter_node2(self):
        global analysis_misfit 
            
        self.reset_vmask()
        self.vmask[4] = True # R2       
        restore_analysis_misfit_function()   
        self.setup()      
        res = run_fitter()             
        r2 = np.abs(res.x[0])           
        print('r2 = '+str(r2))
              
        if res.success:      
            self.R2 = r2           
            self.setup()    # сохраняем вектор с успешным подбором
            
        print(self)
        analysis_plot(title='подбор R2')
    
        
    # запустить подбор значений R3
    def run_fitter_node3(self):
        global analysis_misfit 
            
        self.reset_vmask()        
        self.vmask[7] = True # R3
        
        analysis_misfit = analysis_misfit_by_negativeVoltage
        self.setup()
        
        res = run_fitter()
        
        r3 = np.abs(res.x[0])              
        print('r3 = '+str(r3))
        
        if res.success:            
            self.R3 = r3
            self.setup()    # сохраняем вектор с успешным подбором
            

    # запустить подбор значений R1,R3
    def run_fitter_node13(self):
        global analysis_misfit 
            
        self.reset_vmask()
        self.vmask[0] = True # R1
        self.vmask[7] = True # R3
        
        analysis_misfit = analysis_misfit_by_sko
        self.setup()
        
        res = run_fitter()
        
        r1 = np.abs(res.x[0])
        r3 = np.abs(res.x[1])
        print('r1 = '+str(r1))
        print('r3 = '+str(r3))
        
        if res.success: 
            self.R1 = r1
            self.R3 = r3
            self.setup()    # сохраняем вектор с успешным подбором
            

        print(self)
        analysis_plot(title='подбор R1,R3')
    
    
    def run_fitter(self):
        self.find_initial_value_and_detect_nodes()
        #return 
    
        # рассматриваем частные случаи                  
        # только ветка 1
        if self.node_1==True and self.node_2==False and self.node_3==False:
            self.run_fitter_node1()
            analysis_plot(title='РЕЗУЛЬТАТ')
            return
        
        
        # только ветка 3
        if self.node_3==True and self.node_1==False and self.node_2==False:
            self.run_fitter_node3()
            analysis_plot(title='РЕЗУЛЬТАТ')           
            return
        
        # только ветка 2, нам повезло
        if self.node_2==True and self.node_1==False and self.node_3==False:
            self.run_fitter_node2()
            analysis_plot(title='РЕЗУЛЬТАТ')   
            return
        
        
        # случай симметричной задачи              
        self.run_fitter_node13()
        analysis_plot(title='R1 R3')
        self.run_fitter_node2()    
        analysis_plot(title='R2')
        self.run_fitter_node123()  
        analysis_plot(title='РЕЗУЛЬТАТ')
        return 
    
        self.run_fitter_node13()
        self.run_fitter_node12()
        self.run_fitter_node23() 
        self.run_fitter_node2()            
        self.run_fitter_node13()
        analysis_plot(title='РЕЗУЛЬТАТ')
                        
        
    def run_me(self):        
        # self.init_target1() # тест ветки 1
        # self.reset_values()    
        # self.run_fitter()
          
        # self.init_target2()       
        # self.reset_values()
        # self.run_fitter()
        # return 
    
        self.init_target3()
        self.reset_values()
        self.run_fitter()
        
        
#############################################################################    
    
def main():
    #my_test.test_all()
    #case1 = test_case_1()
    #print(case1)
    Sch = CGeneralCircuit()
    Sch.run_me()
    
    
if __name__=='__main__':
    main()
      
##############################################################################
  
