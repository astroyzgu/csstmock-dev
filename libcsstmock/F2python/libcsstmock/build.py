import cffi  
ffibuilder = cffi.FFI() # 声明了外部函数接口(FFI)对象

with open('mymodule.py') as file_obj:
     content = file_obj.read()
module = content

with open('plugin.h') as file_obj:
     header  = file_obj.read()

ffibuilder.embedding_api(header)
ffibuilder.set_source("my_plugin", r'''
    #include "plugin.h"
''')

ffibuilder.embedding_init_code(module)
#ffibuilder.compile(target= "libplugin.dylib", verbose=True)
#ffibuilder.compile(target= "libcsstplugin.so", verbose=True)
ffibuilder.compile(target= "libcsstplugin.a", verbose=True)
#ffibuilder.compile(target= "libplugin.o", verbose=True)
