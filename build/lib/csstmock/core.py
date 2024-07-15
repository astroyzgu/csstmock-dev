import time

def hello_world(): 
    strftime = time.strftime('%Y%m%d %H:%M:%S', time.localtime(time.time()))
    print( "%s: hello_world"% strftime )
    return "hello_world"

