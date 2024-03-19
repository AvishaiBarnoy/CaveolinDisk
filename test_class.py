class test:
    other_lst = []
    def __init__(self, lst=[]):
        def small(x):
            return [i**2 for i in range(x)]
        self.lst = small(5)
        self.len_lst = len(self.lst)
    def fill_lst(self):
        self.lst.extend(['fill_lst'])


a = test()
print(a, a.lst, a.len_lst)
a.fill_lst()
print(a, a.lst, a.len_lst)

print("lst", a.lst)
print("lst_len", a.len_lst)
