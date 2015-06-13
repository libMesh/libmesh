#!/usr/bin/env python
# Copyright (C) 2015 Akselos

JUST RUN ../../../libmesh_install/capnp/bin/capnp compile -oc++ rb_data.capnp

import os
import subprocess as sp
import shutil

REAL_FILE = 'rb_data_real.capnp'
COMPLEX_FILE = 'rb_data_complex.capnp'
COMPILE_COMMAND = "../../../libmesh_install/capnp/bin/capnp compile -oc++ {}"

if __name__ == "__main__":
    
    with open('rb_data.capnp', 'r') as f:
        schema = f.read()
    
    schema_real = schema.replace("<NUMBER>", "Real")
    schema_complex = schema.replace("<NUMBER>", "Complex")
    
    with open(REAL_FILE, 'w') as f:
        f.write(schema_real)
        
    with open(COMPLEX_FILE, 'w') as f:
        f.write(schema_complex)
        
    sp.check_call(COMPILE_COMMAND.format(REAL_FILE), shell=True)
    sp.check_call(COMPILE_COMMAND.format(COMPLEX_FILE), shell=True)
    
    with open(REAL_FILE + '.c++','r+') as f:
        content = f.read()
        f.seek(0,0)
        f.write("#ifndef LIBMESH_USE_COMPLEX_NUMBERS" + '\n' + content + "#endif") 
    
    with open(COMPLEX_FILE + '.c++','r+') as f:
        content = f.read()
        f.seek(0,0)
        f.write("#ifdef LIBMESH_USE_COMPLEX_NUMBERS" + '\n' + content + "#endif") 
        
    shutil.move(REAL_FILE + '.c++', os.path.join('..', '..', 'src', 'reduced_basis', REAL_FILE + '.C'))
    shutil.move(COMPLEX_FILE + '.c++', os.path.join('..', '..', 'src', 'reduced_basis', COMPLEX_FILE + '.C'))
