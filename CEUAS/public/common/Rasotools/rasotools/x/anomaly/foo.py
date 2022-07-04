import bar

#pythran export foo()
def foo():
    print "foo"
    bar.bar()

foo()
