
# coding: utf-8

# In[74]:

#Loading some libraries
import sys
import numpy as np
sys.path.append("src/data")

from replication.PMotion import Polymer

from replication.tools import load_3D_simus


# In[59]:

def test_simple():
    P = Polymer(0,start=0,end=2,origins=[0])
    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)
    for t in range(3):
        P.increment_time(dt=1,fork_speed=1)

    assert( np.all(P.get_replication_profile(1) ==[1,2,3]))


# In[58]:

def test_simple2():

    P = Polymer(0,start=0,end=2,origins=[0,2])
    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)
    P.add_fork(ptags=["Left1","Right1"],otag=2,new_btags=[None,None],diff_diff_tag=None)

    for t in range(3):
        P.increment_time(dt=1,fork_speed=1)
        
    assert( np.all(P.get_replication_profile(1) ==[1,2,1]))


# In[164]:

def test_symetrie():
    
    eps=1e-4
    P = Polymer(0,start=0,end=2,origins=[0,2])
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
        if t == 1:
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


# In[173]:

def test_remove_passivated1():
    P = Polymer(0,start=0,end=2,origins=[0,2])
    P.add_fork(ptags=["Left","Right"],otag=0,new_btags=[None,None],diff_diff_tag=None)

    for t in range(4):
        #print("Avant",[[p,p.position] for p in P.modules])

        inc = P.increment_time(dt=0.5,fork_speed=1)
        bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = inc
        #print(P.ended[-1].position)
        #print("Apres",[[p,p.position] for p in P.modules])
        #print(P.end)
        if t == 1:
            assert("Left" in alone)
            #print("Here")
            
            for p in P.ended:
                if p.tag in ["Left" ]:
                    #print("path", p.path,len)
                    assert(len(p.path) == 1)
                    
        if t == 4:
            assert(2 in passivated_origin)
   
    assert(np.all(P.get_replication_profile(1) ==[1,2,3]))


# In[178]:

def test_remove_passivated2():
    P = Polymer(0,start=0,end=2,origins=[0,2])
    P.add_fork(ptags=["Left","Right"],otag=2,new_btags=[None,None],diff_diff_tag=None)

    for t in range(4):
        #print("Avant",[[p,p.position] for p in P.modules])

        inc = P.increment_time(dt=0.5,fork_speed=1)
        bind_diff, diff_diff, update_bond, passivated_origin, to_release, alone = inc
        #print(P.ended[-1].position)
        #print("Apres",[[p,p.position] for p in P.modules])
        #print(P.end)
        if t == 1:
            assert("Right" in alone)
            #print("Here")
            
            for p in P.ended:
                if p.tag in ["Right" ]:
                    #print("path", p.path,len)
                    assert(len(p.path) == 1)
                    
        if t == 4:
            assert(2 in passivated_origin)
   
    assert(np.all(P.get_replication_profile(1) ==[3,2,1]))


# In[174]:

test_simple()
test_simple2()
test_symetrie()
test_remove_passivated1()
test_remove_passivated2()

