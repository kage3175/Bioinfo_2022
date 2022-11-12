class test:
    def __init__(self):
        self.n=0
        
def testing(x,class_):
    dict={}
    dict[x]= class_()
    return dict

dict=testing(1,test)

print(dict[1].n)