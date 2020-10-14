def read():
    with open('data.txt') as file:
        tmp = {}
        for line in file:
            line = line.split()
            tmp[line[0]] = line[1]
    return tmp
