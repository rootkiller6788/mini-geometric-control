
import os
BASE = r"F:/nano-everything/mini-complex-control-theory/21. mini-geometric-control/mini-feedback-linearization-geo"

def w(fname, content):
    full = os.path.join(BASE, fname)
    with open(full, 'w') as f:
        f.write(content)
    print(f"  wrote {fname}")

def a(fname, content):
    full = os.path.join(BASE, fname)
    with open(full, 'a') as f:
        f.write(content)
    print(f"  appended to {fname}")
