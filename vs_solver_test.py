# vs_solver_test.py
# номер версии не присвоен
# язык Python
#
# составная часть проекта vs_circuit_solver
# тесты решателя для различный условий
#
# автор В.Симонов, 18-мая-2020
# vasily_simonov@mail.ru, github.com/vasily84
#
# license : это модуль в любом виде можно использовать в любых целях. 
# Ссылка на автора приветствуется, но не является обязательной 
#

from vs_circuit_solver import set_circuit_template,\
    init_target_by_circuitFile, init_target_from_csvFile, \
    generate_circuitFile_by_values,run_fitter,analysis_plot,\
    set_circuit_nominals_and_mask
    
import vs_solver_settings as G

# значения R1,R2
test_A_R1R2 = [1e2, 1e3]


def test_all():
    test_A_all()
    
# выполнить тесты для схемы test_A.cir_t
def test_A_all():
    # секретный метод индийских ученых для запуска тестов
    set_circuit_template('test_A.cir_t')
    set_circuit_nominals_and_mask([1e2,1e2], [True,True] )
    test_A1() 
    G.restore_INIT_param()
    #return 

    test_A2()
    G.restore_INIT_param()
    test_A3()
    G.restore_INIT_param()
    test_A4()
    G.restore_INIT_param()
    test_A5()
    G.restore_INIT_param()
    test_A6()
    G.restore_INIT_param()
    test_A7()
    G.restore_INIT_param()
    test_A8()
    G.restore_INIT_param()
    test_A9()
    G.restore_INIT_param()
    test_A10()
    
    
# тест - автоподбор кривой варьированием R1,R2 для Rcs = 4700
def test_A1():
    print('begin of test A1')
    G.INIT_Rcs = 4700
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_from_csvFile('test_A1_result.csv')
   
    #init_target_by_circuitFile()
    
    if run_fitter(result_cir_file_name='test_A1_result.cir',result_csv_file_name='test_A1_result.csv'):
        analysis_plot(title='test_A1 Rcs=4700',pngName='test_A1_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A1')
    

# тест - автоподбор кривой варьированием R1,R2 для Rcs = 470
def test_A2():
    print('begin of test A2')
    G.INIT_Rcs = 470
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_by_circuitFile()
    
    if run_fitter(result_cir_file_name='test_A2_result.cir',result_csv_file_name='test_A2_result.csv'):
        analysis_plot(title='test_A2 Rcs=470',pngName='test_A2_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A2')
    

# тест - автоподбор кривой варьированием R1,R2 для Rcs = 47
def test_A3():
    print('begin of test A3')
    G.INIT_Rcs = 47
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_by_circuitFile()
    
    if run_fitter(result_cir_file_name='test_A3_result.cir',result_csv_file_name='test_A3_result.csv'):
        analysis_plot(title='test_A3 Rcs=47',pngName='test_A3_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A3')
   

# тест - автоподбор кривой варьированием R1,R2 для разных SNR
def test_A4():
    print('begin of test A4')
    G.INIT_SNR = 20 # моделируем более шумное измерений
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_by_circuitFile()
    
    G.INIT_SNR = 40 # подбираем менее шумный сигнал
    if run_fitter(result_cir_file_name='test_A4_result.cir',result_csv_file_name='test_A4_result.csv'):
        analysis_plot(title='test_A4 SNR1 = 20, SNR2 = 40',pngName='test_A4_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A4')
    

# тест - автоподбор кривой варьированием R1,R2 для разных SNR
def test_A5():
    print('begin of test A5')
    G.INIT_SNR = 20 # моделируем более шумное измерений
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_by_circuitFile()
    
    G.INIT_SNR = 20 # подбираем менее шумный сигнал
    if run_fitter(result_cir_file_name='test_A5_result.cir',result_csv_file_name='test_A5_result.csv'):
        analysis_plot(title='test_A5 SNR1 = 20, SNR2 = 20',pngName='test_A5_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A5')
    

# тест - подбор кривой с большим числом точек
def test_A6():
    if G.USE_LIBIVCMP:
        print('test_A6 impossible for libivcmp lib')
        return 
    
    print('begin of test A6')
    G.MAX_NUM_POINTS = 10000
    
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_by_circuitFile()
    
    if run_fitter(result_cir_file_name='test_A6_result.cir',result_csv_file_name='test_A6_result.csv'):
        analysis_plot(title='test_A6 MAX_NUM_POINTS = 10000',pngName='test_A6_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A6')


# тест - подбор кривой c малым числом точек 
def test_A7():
    if G.USE_LIBIVCMP:
        print('test_A7 impossible for libivcmp lib')
        return 
    
    print('begin of test A7')
    G.MAX_NUM_POINTS = 10
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_by_circuitFile()
    
    if run_fitter(result_cir_file_name='test_A7_result.cir',result_csv_file_name='test_A7_result.csv'):
        analysis_plot(title='test_A7 MAX_NUM_POINTS = 10',pngName='test_A7_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A1')
    
    
# тест - подбор кривой с большей амплитудой напряжения
def test_A8():
    print('begin of test A8')  
    G.INIT_V = 3
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_by_circuitFile()
    G.INIT_V = 3.3
    if run_fitter(result_cir_file_name='test_A8_result.cir',result_csv_file_name='test_A8_result.csv'):
        analysis_plot(title='test_A8 V1=4 Volt V2 = 3.3 Volt',pngName='test_A8_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A8')
    
    
# тест - подбор кривой с большей значением Rcs 
def test_A9():
    print('begin of test A9')
    G.INIT_Rcs = 4700
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_by_circuitFile()
    
    G.INIT_Rcs = 3300
    if run_fitter(result_cir_file_name='test_A9_result.cir',result_csv_file_name='test_A9_result.csv'):
        analysis_plot(title='test_A9 Rcs1=4700 Ohm Rcs2 = 3300 Ohm',pngName='test_A9_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A9')
    
    
# тест - подбор кривой с удвоенным значением Rcs и напряжением
def test_A10():
    print('begin of test A10')
    G.INIT_Rcs = 4700
    G.INIT_V = 3
    generate_circuitFile_by_values(test_A_R1R2)
    init_target_by_circuitFile()
    
    G.INIT_Rcs = 2.*4700
    G.INIT_V = 6
    if run_fitter(result_cir_file_name='test_A10_result.cir',result_csv_file_name='test_A10_result.csv'):
        analysis_plot(title='test_A10 Rcs and Voltage multiplied by 2',pngName='test_A10_result.png')
    else:
        print('!!! OPTIMIZATION FAILED !!!')
        
    print('end of test A10')
