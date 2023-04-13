#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sy             #importation de sympy
sy.init_printing()             #affichage en latex des outputs
x=sy.Symbol('x')
def Berlekamp(seq):
    f=1
    deg_f_list=[]
    n=len(seq)
    f_list=[1]
    for i in range(n):
        #print("Step",i+1)
        d=0
        deg_f=sy.degree(f)
        deg_f_list.append(deg_f)
        coeffs_f=sy.Poly(f,x).all_coeffs()
        for j in range(deg_f+1):
            d=d+coeffs_f[deg_f-j] * seq[i-j]
        if d==0:
            f=sy.Poly(f,x).trunc(2)
            f_list.append(f)
        else:
            L_index=[k for k in range(len(deg_f_list)-1) if deg_f_list[k]<deg_f_list[k+1]]
            if len(L_index)==0:
                m=-1
                fm=1
                Lm=0
            else:
                m=max(L_index)
                fm=f_list[m]
                Lm=deg_f_list[m]
            #print("deg of f = ", deg_f)
            #print("m = ",m)
            #print("list of f = ",f_list)
                Li=deg_f_list[i]
                if (m-Lm)>=(i-Li):
                    f=f+x**((m-Lm)-(i-Li))*fm
                    f=sy.Poly(f,x).trunc(2)
                    f_list.append(f)
                else:
                    f=(x**((i-Li)-(m-Lm)))*(f)+fm
                    f=sy.Poly(f,x).trunc(2)
                    f_list.append(f)
        #print(deg_f_list)    
        #print(f)
    return f.trunc(2)


# In[2]:


seq=[1, 1, 0, 1, 0, 0, 1]
Berlekamp(seq)


# In[ ]:




