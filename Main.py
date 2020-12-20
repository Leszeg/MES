import MES
import time

if __name__ == '__main__':
    start = time.time()
    MES.run()
    end = time.time()
    print(f"Czas wykonania {end - start}")
