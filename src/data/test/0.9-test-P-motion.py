
# coding: utf-8

# In[1]:

#Loading some libraries

from pylab import *
import numpy as np

import copy
import sys
sys.path.append("src/data")
from replication.PMotion import Polymer
from replication.simulate_1D import simulate

from replication.tools import load_3D_simus


# In[2]:

def test_simple():
    P = Polymer(0,start=0,end=2,origins=[0])
    #
    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)

    for t in range(3):
        #print([m.path for m in P.modules if m.move])

        P.increment_time(dt=1,fork_speed=1)
    #print(P.get_replication_profile(1))
    assert( np.all(P.get_replication_profile(1) ==[1,2,3]))


# In[3]:

def test_simple2():

    P = Polymer(0,start=0,end=2,origins=[0,2])
    P.increment_time(dt=0.5,fork_speed=1)
    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)
    P.add_fork(ptags=["Left1","Right1"],otag=2,new_btags=[None,None],diff_diff_tag=None)

    for t in range(3):
        P.increment_time(dt=1,fork_speed=1)
    
    #print([m.path for m in P.ended if m.move])
    #print(P.get_replication_profile(1))
    assert( np.all(P.get_replication_profile(1) ==[1,2,1]))
#test_simple2()


# In[ ]:




# In[4]:

def test_symetrie():
    
    eps=1e-4
    P = Polymer(0,start=0,end=2,origins=[0,2])
    P.increment_time(dt=0.5,fork_speed=1)

    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)
    P.add_fork(ptags=["Left1","Right1"],otag=2,new_btags=[None,None],diff_diff_tag=None)

    for t in range(4):
        #print("Avant",[[p,p.position] for p in P.modules])

        inc = P.increment_time(dt=0.5,fork_speed=1)
        #print(inc)
        bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = inc
        #print(P.ended[-1].position)
        #print("Apres",[[p,p.position] for p in P.modules])
        #print(P.end)
        if t == 0:
            assert("Right1" in alone and "Left" in alone)
            #print("Here")
            for p in P.ended:
                if p.tag in ["Left" ,"Right1"]:
                    #print("path", p.path,len)
                    assert(len(p.path) == 1)
                elif p.tag in ["Left1" ,"Right"]:
                    assert(len(p.path) == 2)
                    if p.tag == "Right":
                        #print(p.path)
                        assert( p.path[0].pos == 0.5 and p.path[1].pos == 1)
                    if p.tag == "Left1":
                        #print(p.path,p.position)
                        assert( p.path[0].pos == 2.5 and p.path[1].pos == 1.5)
                    #print(p.tag,p.path)

    assert(np.all(P.get_replication_profile(1) ==[1,2,1]))


# In[5]:

def test_remove_passivated1():
    P = Polymer(0,start=0,end=2,origins=[0,2])
    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)

    for t in range(5):
        #print("Avant",[[p,p.position] for p in P.modules])

        inc = P.increment_time(dt=0.5,fork_speed=1)
        bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = inc
        #print(P.ended[-1].position)
        #print("Apres",[[p,p.position] for p in P.modules])
        #print(P.end)
        if t == 0:
            assert("Left" in alone)
            #print("Here")
            
            for p in P.ended:
                if p.tag in ["Left" ]:
                    #print("path", p.path,len)
                    assert(len(p.path) == 1)
                    
        if t == 4:
            assert(2 in passivated_origin)
   
    assert(np.all(P.get_replication_profile(1) ==[1,2,3]))
test_remove_passivated1()


# In[6]:

def test_remove_passivated2():
    P = Polymer(0,start=0,end=2,origins=[0,2])
    P.add_fork(ptags=["Left","Right"],otag=2,new_btags=[None,None],diff_diff_tag=None)

    for t in range(5):
        #print("Avant",[[p,p.position] for p in P.modules])

        inc = P.increment_time(dt=0.5,fork_speed=1)
        bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = inc
        #print(P.ended[-1].position)
        #print("Apres",[[p,p.position] for p in P.modules])
        #print(P.end)
        #print(t,passivated_origin)
        if t == 0:
            assert("Right" in alone)
            #print("Here")
            
            for p in P.ended:
                if p.tag in ["Right" ]:
                    #print("path", p.path,len)
                    assert(len(p.path) == 1)
                    
        if t == 4:
            assert(0 in passivated_origin)
    assert(np.all(P.get_replication_profile(1) ==[3,2,1]))
#test_remove_passivated2()


# In[7]:

def test_remove_passivated3():
    P = Polymer(0,start=0,end=2,origins=[0,1,2])
    P.modules[0].position = 0.9
    P.modules[1].position = 1.1
    
    P.add_fork(ptags=["Left","Right"],otag=2,new_btags=[None,None],diff_diff_tag=None)

    for t in range(5):
        #print("Avant",[[p,p.position] for p in P.modules])

        inc = P.increment_time(dt=1,fork_speed=1)
        bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = inc
        #print(P.ended[-1].position)
        #print("Apres",[[p,p.position] for p in P.modules])
        #print(P.end)
        #print(P.modules)

        #print(t,passivated_origin)
        if t == 0:
            assert("Right" in alone)
            #print("Here")
            
            for p in P.ended:
                if p.tag in ["Right" ]:
                    #print("path", p.path,len)
                    assert(len(p.path) == 1)
                    
        if t == 1:
            assert(0 in passivated_origin and  1 in passivated_origin)
    #print(P.get_replication_profile(1))
    assert(np.all(P.get_replication_profile(1) == [3,2,1]))
#test_remove_passivated3()


# In[8]:

def test_remove_passivated4():
    P = Polymer(0,start=0,end=2,origins=[0,1,2])
    P.modules[1].position = 1.9
    P.modules[2].position = 2.1
    
    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)

    for t in range(5):
        #print("Avant",[[p,p.position] for p in P.modules])

        inc = P.increment_time(dt=1,fork_speed=1)
        bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = inc
        #print(P.ended[-1].position)
        #print("Apres",[[p,p.position] for p in P.modules])
        #print(P.end)
        #print(P.modules)
        #print(t,passivated_origin)
        if t == 0:
            assert("Left" in alone)
            #print("Here")
            
            for p in P.ended:
                if p.tag in ["Right" ]:
                    #print("path", p.path,len)
                    assert(len(p.path) == 1)
                    
        if t == 1:
            assert(1 in passivated_origin and  2 in passivated_origin)
    #print(P.get_replication_profile(1))
    assert(np.all(P.get_replication_profile(1) ==[1,2,3]))
#test_remove_passivated2()
#test_remove_passivated4()


# In[9]:

def test_remove_passivated5():
    P = Polymer(0,start=0,end=2,origins=[0,1,2])
    P.modules[1].position = 1.9
    P.modules[2].position = 2.1
    
    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)

    for t in range(5):
        #print("Avant",[[p,p.position] for p in P.modules])

        inc = P.increment_time(dt=0.5,fork_speed=2.)
        bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = inc
        #print(P.ended[-1].position)
        #print("Apres",[[p,p.position] for p in P.modules])
        #print(P.end)
        #print(P.modules)
        #print(t,passivated_origin)
        if t == 0:
            assert("Left" in alone)
            #print("Here")
            
            for p in P.ended:
                if p.tag in ["Right" ]:
                    #print("path", p.path,len)
                    assert(len(p.path) == 1)
                    
        if t == 1:
            assert(1 in passivated_origin and  2 in passivated_origin)
            
    #print(P.get_replication_profile(2))
    #print(P.g)
    assert(np.all(P.get_replication_profile(2) ==[1,1,2]))
#test_remove_passivated2()
#test_remove_passivated5()


# In[10]:

def test_jump2():
    P = Polymer(0,start=0,end=2,origins=[0,2])
    #
    P.add_fork(ptags=["Left","Right"],otag=2,new_btags=[None,None],diff_diff_tag=None)

    for t in range(3):
        #print([m.path for m in P.modules if m.move])

        P.increment_time(dt=1,fork_speed=2.6)
    assert( np.all(P.get_replication_profile(2.6) ==[1,1,1]))


# In[11]:

def test_jump():
    P = Polymer(0,start=0,end=2,origins=[0])
    #
    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)

    for t in range(3):
        #print([m.path for m in P.modules if m.move])

        P.increment_time(dt=1,fork_speed=2.6)
    #print(P.get_DNA_with_time())
    assert( np.all(P.get_replication_profile(2.6) ==[1,1,1]))


# In[12]:

def test_fork_speed_dt_speed():
    testv = [[1,1],[1,0.5],[1,2],[2,1],[2,0.5],[0.5,1],[0.5,0.5]]
    for fork_speed,dt_speed,in testv:

        #def test_simulate_1d():
        for random in [True,False]:
            S = simulate(nori=[np.arange(20)],ndiff=2,lengths=[20],p_on=0.5,
                         p_off=0,only_one=True,fork_speed=fork_speed,dt_speed=dt_speed,gindin=False,random=random)
            S.simulate(400)

            assert(S.polys[0].get_DNA_with_time(fork_speed)[0][-1] == 20)
            assert(S.polys[0].get_DNA_with_time(fork_speed,dt=dt_speed)[0][-1] == 20)

            assert(S.polys[0].get_DNA_with_time(fork_speed)[0][0] == 0)


# In[13]:

test_simple()
test_simple2()
test_symetrie()
test_remove_passivated1()
test_remove_passivated2()
test_remove_passivated3()
test_remove_passivated4()
test_jump()
test_jump2()
test_fork_speed_dt_speed()


# In[ ]:



