# vs_circuit_solver.py
# номер версии не присвоен
# язык Python
#
# модуль подбора значений R,C для вариантов электронной схемы
# исходя из моделирования подобной схемы в ngspice
# поставляется без всякой оптимизации, ибо имеет целью установление методики 
# расчета таких вещей и определения границ применимости этой методики
#
# автор В.Симонов, 14-мая-2020
# vasily_simonov@mail.ru, github.com/vasily84
#
# license : это модуль в любом виде можно использовать в любых целях. 
# Ссылка на автора приветствуется, но не является обязательной 
#

import scipy.optimize as spo
import numpy as np
import matplotlib.pyplot as plt
from ctypes import c_double

# модуль Михаила Лукьянова, github.com/LukyanovM/MySpice 
import MySpice as spice

# закомментировать импорт libivcmp, если этой библиотеки нет
import libivcmp

# использовать внешнюю библиотеку libivcmp для сравнения кривых?
# !! установить True или False 
USE_LIBIVCMP = False

# метод оптимизации функции подбора параметров R,C,
# варианты для функции scipy.optimize.minimum() 
# !! Раскомментировать необходимый FITTER_METHOD
FITTER_METHOD = 'Powell' # это метод работает лучше всего
#FITTER_METHOD = 'Nelder-Mead' # тоже рабочий
#FITTER_METHOD = None # не рабочий - происходит потеря точности 

# погрешность подбора кривых
TOLERANCE = 1e-7
# число точек в массивах тока и напряжения
MAX_NUM_POINTS = 1000

if USE_LIBIVCMP:
    TOLERANCE = 1e-15
    MAX_NUM_POINTS = libivcmp.MAX_NUM_POINTS


# частота, Гц
TARGET_F = 1000
# амплитудное напряжение, В
TARGET_V = 3
# токоограничивающий резистор, Ом
TARGET_Rcs = 4700
# SIGNAL/NOISE ratio
TARGET_SNR = 40.0

    
# глобальные переменные для работы PySpice
circuit = None
input_data = None
analysis = None

# целевая кривая с током. Та, которую мы подбираем.
target_VCurrent = None
target_input_dummy = None

# целевая кривая с током для сравнения в библиотеке libivcmp
target_IVCurve = None

#глобальные переменные
circuitTemplateFileName = 'test_A.cir_t'
circuitSessionFileName = 'var1.cir'


# инициализировать целевую модель
def init_target_by_circuitFile(fileName):
    global target_VCurrent, target_input_dummy, target_IVCurve
    process_circuitFile(fileName)
    target_VCurrent = analysis.VCurrent
    target_input_dummy = analysis.input_dummy
    
    if USE_LIBIVCMP:
        target_IVCurve = analysis_to_IVCurve()
        
# последний анализ перевести в форму, пригодную для сравниния в libivcmp
def analysis_to_IVCurve():
    iv_curve = libivcmp.IvCurve()
    for i in range(MAX_NUM_POINTS):
        iv_curve.voltages[i] = c_double(analysis.VCurrent[i])
        iv_curve.currents[i] = c_double(analysis.input_dummy[i])
        
    libivcmp.SetMinVC(0, 0)
    return iv_curve
    
                
# промоделировать файл схемы
def process_circuitFile(fileName='default1.cir',csvName=''):
    global circuit,input_data,analysis
    circuit = spice.LoadFile(fileName)
    input_data = spice.Init_Data(TARGET_F, TARGET_V, TARGET_Rcs,TARGET_SNR )
    analysis = spice.CreateCVC1(circuit, input_data, MAX_NUM_POINTS, "input", 10)
    if(not csvName==''):
        spice.SaveFile(analysis, csvName)
 
    
# вывести на график результат моделирования
def analysis_plot(title='',pngName=''):
    global circuit,input_data,analysis
    figure1 = plt.figure(1, (20, 10))
    plt.grid()
    plt.plot(target_input_dummy, target_VCurrent,color='red')
    plt.plot(analysis.input_dummy, analysis.VCurrent,color='blue')
    if (not title==''):
       plt.title(title)       
    plt.xlabel('Напряжение [В]')
    plt.ylabel('Сила тока [А]')
    if(not pngName==''):
        plt.savefig(pngName)
        
    plt.show()
    
    
       
# считать файл шаблона схемы, сделать замену {} на значения варьирования Xi_values,
# сохранить с новым именем
def generate_circuitFile_by_values(fileTemplate, fileCircuit, Xi_values):
    templateF = open(fileTemplate)
    newF = open(fileCircuit, 'w')
    i = 0
    for tStr in templateF:
        cStr = tStr
        if tStr.find('{}')>=0:
            cStr = tStr.format(str(Xi_values[i]))
            i += 1
        newF.write(cStr)
                
    templateF.close()
    newF.close()


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

    
if USE_LIBIVCMP:
    analysis_misfit = analysis_misfit_by_libivcmp
else:
    analysis_misfit = analysis_misfit_by_sko
    
    
# для варьирования без ограничений, специально выбранная функция
def my_abs(v):
    if v>=0.0:
        return v 
    return 2.*abs(v)+v*v


# счетчик числа вызова функции оптимизатором
ffCount = 0
# функция вызывается оптимизатором
def fitter_subroutine(Xargs):
    global ffCount
    x1 = np.abs(Xargs)
    generate_circuitFile_by_values(circuitTemplateFileName,circuitSessionFileName,x1)
    process_circuitFile(fileName=circuitSessionFileName)
    ffCount += 1
    print('fitter_subroutine Count = '+str(ffCount))
    print(x1)
    m = analysis_misfit()
    print('misfit = '+str(m))
    return m


# запустить автоподбор
def run_fitter(result_cir_file_name='',result_csv_file_name=''):
    global ffCount
    ffCount = 0
    maxV = np.amax(target_VCurrent)
    maxI = np.amax(target_input_dummy)
    r = abs(maxV)/(abs(maxI)+0.01) # стартовые значения для варьирования  
    Xargs = [2.*r,2.*r]
    
    resX = spo.minimize(fitter_subroutine,Xargs,method=FITTER_METHOD,tol=TOLERANCE,options={'maxiter':100})
    # вызываем с результатом оптимизации, ибо предыдущий вызов может быть неоптимальным
    fitter_subroutine(resX.x) 
    print(resX.message)
    print('result X = '+str(resX.x))
    print('function evaluation Count = '+str(ffCount))
    print('fitter method N iteration = '+str(resX.nit))
    
    if resX.success:
        if(not result_csv_file_name==''):
            spice.SaveFile(analysis, result_csv_file_name)
        if(not result_cir_file_name==''):
            x1 = np.abs(resX.x)
            generate_circuitFile_by_values(circuitTemplateFileName,result_cir_file_name,x1)
    
    return resX.success

##############################################################################
# серия тестов для схемы test_A.cir_t
    
_TOLERANCE = TOLERANCE
_MAX_NUM_POINTS = MAX_NUM_POINTS 
_TARGET_F = TARGET_F
_TARGET_V = TARGET_V
_TARGET_Rcs = TARGET_Rcs
_TARGET_SNR = TARGET_SNR

def restore_TARGET_param():
    global TOLERANCE,MAX_NUM_POINTS,TARGET_F,TARGET_V,TARGET_Rcs,TARGET_SNR
    TOLERANCE = _TOLERANCE
    MAX_NUM_POINTS = _MAX_NUM_POINTS 
    TARGET_F = _TARGET_F
    TARGET_V = _TARGET_V
    TARGET_Rcs = _TARGET_Rcs
    TARGET_SNR = _TARGET_SNR
    

# выполнить тесты для схемы test_A.cir_t
def test_A_all():
    #return 
    # секретный метод индийских ученых для запуска тестов
    test_A1() 
    restore_TARGET_param()
    test_A2()
    restore_TARGET_param()
    test_A3()
    restore_TARGET_param()
    test_A4()
    restore_TARGET_param()
    test_A5()
    restore_TARGET_param()
    test_A6()
    restore_TARGET_param()
    test_A7()
    restore_TARGET_param()
    test_A8()
    restore_TARGET_param()
    test_A9()
    restore_TARGET_param()
    test_A10()


# тест - автоподбор кривой варьированием R1,R2 для Rcs = 4700
def test_A1():
    global TARGET_Rcs,circuitTemplateFileName
    print('begin of test A1')
    TARGET_Rcs = 4700
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    
    if run_fitter(result_cir_file_name='test_A1_result.cir',result_csv_file_name='test_A1_result.csv'):
        analysis_plot(title='test_A1 Rcs=4700',pngName='test_A1_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A1')
    

# тест - автоподбор кривой варьированием R1,R2 для Rcs = 470
def test_A2():
    global TARGET_Rcs,circuitTemplateFileName
    print('begin of test A2')
    TARGET_Rcs = 470
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    
    if run_fitter(result_cir_file_name='test_A2_result.cir',result_csv_file_name='test_A2_result.csv'):
        analysis_plot(title='test_A2 Rcs=470',pngName='test_A2_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A2')
    

# тест - автоподбор кривой варьированием R1,R2 для Rcs = 47
def test_A3():
    global TARGET_Rcs,circuitTemplateFileName
    print('begin of test A3')
    TARGET_Rcs = 47
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    
    if run_fitter(result_cir_file_name='test_A3_result.cir',result_csv_file_name='test_A3_result.csv'):
        analysis_plot(title='test_A3 Rcs=47',pngName='test_A3_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A3')
   

# тест - автоподбор кривой варьированием R1,R2 для разных SNR
def test_A4():
    global TARGET_SNR,circuitTemplateFileName
    print('begin of test A4')
    TARGET_SNR = 20 # моделируем более шумное измерений
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    
    TARGET_SNR = 40 # подбираем менее шумный сигнал
    if run_fitter(result_cir_file_name='test_A4_result.cir',result_csv_file_name='test_A4_result.csv'):
        analysis_plot(title='test_A4 SNR1 = 20, SNR2 = 40',pngName='test_A4_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A4')
    

# тест - автоподбор кривой варьированием R1,R2 для разных SNR
def test_A5():
    global TARGET_SNR,circuitTemplateFileName
    print('begin of test A5')
    TARGET_SNR = 20 # моделируем более шумное измерений
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    
    TARGET_SNR = 20 # подбираем менее шумный сигнал
    if run_fitter(result_cir_file_name='test_A5_result.cir',result_csv_file_name='test_A5_result.csv'):
        analysis_plot(title='test_A5 SNR1 = 20, SNR2 = 20',pngName='test_A5_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A5')
    

# тест - подбор кривой с большим числом точек
def test_A6():
    global circuitTemplateFileName,MAX_NUM_POINTS
    if USE_LIBIVCMP:
        print('test_A6 impossible for libivcmp lib')
        return 
    
    print('begin of test A6')
    MAX_NUM_POINTS = 10000
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    
    if run_fitter(result_cir_file_name='test_A6_result.cir',result_csv_file_name='test_A6_result.csv'):
        analysis_plot(title='test_A6 MAX_NUM_POINTS = 10000',pngName='test_A6_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A6')


# тест - подбор кривой c малым числом точек 
def test_A7():
    global MAX_NUM_POINTS,circuitTemplateFileName
    if USE_LIBIVCMP:
        print('test_A7 impossible for libivcmp lib')
        return 
    
    print('begin of test A7')
    MAX_NUM_POINTS = 10
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    
    if run_fitter(result_cir_file_name='test_A7_result.cir',result_csv_file_name='test_A7_result.csv'):
        analysis_plot(title='test_A7 MAX_NUM_POINTS = 10',pngName='test_A7_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A1')
    
    
# тест - подбор кривой с большей амплитудой напряжения
def test_A8():
    global circuitTemplateFileName,TARGET_V,MAX_NUM_POINTS
    print('begin of test A8')  
    TARGET_V = 3
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    TARGET_V = 3.3
    if run_fitter(result_cir_file_name='test_A8_result.cir',result_csv_file_name='test_A8_result.csv'):
        analysis_plot(title='test_A8 V1=4 Volt V2 = 3.3 Volt',pngName='test_A8_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A8')
    
    
# тест - подбор кривой с большей значением Rcs 
def test_A9():
    global TARGET_Rcs,circuitTemplateFileName
    print('begin of test A9')
    TARGET_Rcs = 4700
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    
    TARGET_Rcs = 3300
    if run_fitter(result_cir_file_name='test_A9_result.cir',result_csv_file_name='test_A9_result.csv'):
        analysis_plot(title='test_A9 Rcs1=4700 Ohm Rcs2 = 3300 Ohm',pngName='test_A9_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A9')
    
    
# тест - подбор кривой с удвоенным значением Rcs и напряжением
def test_A10():
    global TARGET_Rcs,circuitTemplateFileName,TARGET_V
    print('begin of test A10')
    TARGET_Rcs = 4700
    TARGET_V = 3
    circuitTemplateFileName = 'test_A.cir_t'
    generate_circuitFile_by_values(circuitTemplateFileName,'var1.cir',[10,100])
    init_target_by_circuitFile('var1.cir')
    
    TARGET_Rcs = 2.*4700
    TARGET_V = 6
    if run_fitter(result_cir_file_name='test_A10_result.cir',result_csv_file_name='test_A10_result.csv'):
        analysis_plot(title='test_A10 Rcs and Voltage multiplied by 2',pngName='test_A10_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A10')
    
       
##############################################################################
def main():
    test_A_all()
    
if __name__=='__main__':
    main()
    
    
##############################################################################
  
