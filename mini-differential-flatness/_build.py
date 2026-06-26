import os, base64
def w(path, b64):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(base64.b64decode(b64).decode())
    print(f"Wrote {path}")
base = "."