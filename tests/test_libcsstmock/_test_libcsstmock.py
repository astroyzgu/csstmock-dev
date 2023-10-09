import os 
def test_libcsstmock():  
    assert 0 == os.system('gfortran test4.f -o test4 -lcsstmock')
    assert 0 == os.system('test4')
    assert 0 == os.system('gfortran test0.f -o test0 -lcsstmock')
    assert 0 == os.system('test0')
