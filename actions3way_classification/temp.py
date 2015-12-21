def run():
    for i in xrange(3):
        for j in xrange(4):
            yield i+1, i+2

a = run()
print "bla"