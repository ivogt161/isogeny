for i in xrange(2, 7):
    N = 2**i
    for line in open("gl_2_full/gl_2_" + str(N) + ".dat", "r"):
        label = line.split()[0]
        spl = line.split(":")
        if len(spl) > 1:
            gens = line.split(":")[1]
            gen = gens[1:-2]
            g = open("gl2_" + str(N) + ".txt", "a")
            g.write(str(label) + ":" + "[" + str(gen) + "]" + "\n")
            g.close()
