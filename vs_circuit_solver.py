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
Xi_long = [0.,0.]

# Маска оптимизируемых параметров - список булевого типа, например -
# Xi_long = [10.0, x1, 300., x3]
# Xi_mask = [False,True,False,True] -> X_short = [x1,x3]
Xi_mask = [True,True]


def Xi_unroll(x_short):
    global Xi_long
    XL = Xi_long
    
    j = 0
    for i in range(0,len(Xi_mask)):
        if Xi_mask[i]:
            XL[i] = x_short[j]
            j += 1
            
    return XL


# установить все известные и номиналы и маску оптимизации
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
            cStr = tStr.format(str(Xi_values[i]))
            i +=1
            
        newF.write(cStr)
    newF.close()

        
# промоделировать файл схемы
def process_circuitFile(fileName=circuit_SessionFileName,csvName=''):
    global analysis
    circuit = spice.LoadFile(fileName)
    input_data = spice.Init_Data(G.INIT_F, G.INIT_V, G.INIT_Rcs,G.INIT_SNR )
    analysis = spice.CreateCVC1(circuit, input_data, G.MAX_NUM_POINTS, "input", 10)
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
    plt.plot(firstStep_input_dummy,firstStep_VCurrent,color='yellow')
    if (not title==''):
       plt.title(title)       
    plt.xlabel('Напряжение [В]')
    plt.ylabel('Сила тока [А]')
    if(not pngName==''):
        plt.savefig(pngName)
        
    plt.show()
                           

# вычислить несовпадение последнего анализа и целевой функции. 
# возвращает неотрицательное число
def analysis_misfit_by_sko():   
    s = target_VCurrent-analysis.VCurrent    
    #s2 = s*s
    s2 = np.abs(s) # работает наилучшим образом
    s_sum = np.sum(s2)
    return s_sum


# вычислить несовпадение последнего анализа и целевой функции. 
# использует модуль libivcmp
def analysis_misfit_by_libivcmp():
    step_IVCurve = analysis_to_IVCurve()
    res = libivcmp.CompareIvc(target_IVCurve, step_IVCurve, libivcmp.MAX_NUM_POINTS)
    return res

# установить используемую функция вычисления несовпадения
if G.USE_LIBIVCMP:
    analysis_misfit = analysis_misfit_by_libivcmp
else:
    analysis_misfit = analysis_misfit_by_sko
    
    
# для варьирования без ограничений, специально выбранная функция
def my_abs(vect):
    vect2 = []
    for v in vect:
        if v>=0.0:
            vect2 +=[v]
        else:
            vect2 +=[2.*abs(v)+v*v]
            
    return vect2


# счетчик числа вызова функции оптимизатором
ffCount = 0
# функция вызывается оптимизатором
def fitter_subroutine(Xargs):
    global ffCount
    x1 = np.abs(Xargs)
    #x1 = my_abs(Xargs)
    generate_circuitFile_by_values(Xi_unroll(x1))
    process_circuitFile()
    ffCount += 1
    print('fitter_subroutine Count = '+str(ffCount))
    print(x1)
    m = analysis_misfit()
    print('misfit = '+str(m))
    return m


# найти стартовые значения для варьирования
def find_init_Xi():
    maxV = np.amax(target_VCurrent)
    maxI = np.amax(target_input_dummy)
    r = abs(maxV)/(abs(maxI)+0.01)   
    
    return [2.*r, 2.*r]


# запустить автоподбор
def run_fitter(result_cir_file_name='',result_csv_file_name=''):
    global ffCount, firstStep_VCurrent, firstStep_input_dummy
    ffCount = 0
    Xargs = find_init_Xi()
    
    fitter_subroutine(Xargs)
    firstStep_VCurrent = analysis.VCurrent
    firstStep_input_dummy = analysis.input_dummy
    
    resX = spo.minimize(fitter_subroutine,Xargs,method=G.FITTER_METHOD,tol=G.TOLERANCE,options={'maxiter':100})
    
    #bounds = [(-1., 1e3), (-1., 1e3)]
    #resX = spo.shgo(fitter_subroutine,bounds)
    #resX = spo.dual_annealing(fitter_subroutine,bounds)
    #resX = spo.differential_evolution(fitter_subroutine,bounds)
    #resX = spo.basinhopping(fitter_subroutine,bounds)
    # вызываем с результатом оптимизации, ибо предыдущий вызов может быть неоптимальным
    fitter_subroutine(resX.x) 
    print('\n')
    print(resX.message)
    print('result X = '+str(resX.x))
    print('function evaluation Count = '+str(ffCount))
    print('fitter method N iteration = '+str(resX.nit))
    
    if resX.success:
        if(not result_csv_file_name==''):
            spice.SaveFile(analysis, result_csv_file_name)
        if(not result_cir_file_name==''):
            x1 = np.abs(resX.x)
            generate_circuitFile_by_values(x1,result_cir_file_name,)
    
    return resX.success
            
       
##############################################################################
    
def main():
    #init_target_from_csvFile('test_A1_result.csv')
    my_test.test_all()
    
if __name__=='__main__':
    main()
      
##############################################################################
  
