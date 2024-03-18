class test:
    other_lst = []
    def __init__(self, lst=[]):
        def small(x):
            return [i**2 for i in range(x)]
        self.lst = small(5)

a = test()
print(a, a.lst)
a.fill_lst()
print(a, a.lst)
