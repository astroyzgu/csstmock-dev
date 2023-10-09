import csstmock
from csstmock import core
print(csstmock.__path__) 

def test_basic(): 
    assert "hello_world" == core.hello_world() 

