# -*- coding: utf-8 -*-
"""
Created on Wed Apr 01 14:55:55 2015

@author: Patrick Kilecdi
This module stores information for the space group
sgroup.wycgen generate a random fractional coordinate sitting on a Wyckoff position
sgroup.wycarr produces a list of Wyckoff position arrangements that add up to nmpc of molecules per cell
"""
import random
from utilities.write_log import print_time_log

setting={0:random.random(),1:0.25,2:0.75}
enantiomorphic=set([3,4,5,6,7,8,9,10,11,12,13,14,15])

MAX_AVAILABLE_SPACE_GROUP = 142
CHIRAL_SPACE_GROUPS = \
    ([1] + range(3,6) + range(16,25) + range(75,81) + range(89,99) 
    + range(143,147) + range(149,156) + range(168,174) + range(177,183)
    + range(195,200) + range(207,215))
RACEMIC_SPACE_GROUPS = \
    [x for x in range(1, 231) if not x in CHIRAL_SPACE_GROUPS] 

class SpaceGroupManager():
    def __init__(
        self, nmpc, is_chiral, wyckoff_list=[0],
        space_groups_allowed=None):
        '''
        nmpc: Number of molecule per cell; allowed space group must have Wyckoff
            position combo that can match nmpc
        is_chiral: Allowed space groups must be either chiral or racemic
        wyckoff_list: Constrain the Wyckoff position selection.
            Default [0], molecules must be placed on a general Wyckoff position.
            Set this to None to allow any wyckoff position combo
        space_groups_allowed: A list of int with allowed space groups; will be
            further pruned to be compatible with the above

        Raises ValueError if no valid space groups are allowed
        '''
        self._maximum_space_group = MAX_AVAILABLE_SPACE_GROUP
        self._chiral_space_groups = CHIRAL_SPACE_GROUPS
        self._racemic_space_groups = RACEMIC_SPACE_GROUPS
        self._nmpc = nmpc
        self._is_chiral = is_chiral
        self._wyckoff_list = wyckoff_list
        self._space_groups_allowed = space_groups_allowed
        self._deduce_allowed_space_groups()

        self._space_group = None
        self._space_group_selected_counter = {}
        self._space_group_user_counter = {}
        for sg in self._space_groups_allowed:
            self._space_group_selected_counter[sg] = 0
            self._space_group_user_counter[sg] = 0

    def _deduce_allowed_space_groups(self):
        space_group_range = []
        if self._space_groups_allowed != None:
            space_group_range = self._space_groups_allowed
        else:
            space_group_range = range(0, self._maximum_space_group + 1)

        if self._is_chiral:
            space_group_range = \
                [sg for sg in space_group_range
                    if sg in self._chiral_space_groups]
        else:
            space_group_range = \
                [sg for sg in space_group_range
                    if sg in self._racemic_space_groups]
        
        self._space_group_range = []
        if (self._wyckoff_list != None):
            self._is_wyckoff_list_fixed = True
            for sg in space_group_range:
                space_group = Sgroup(sg)
                if (space_group.wyckoff_counter(self._wyckoff_list)
                    == self._nmpc):
                    self._space_group_range.append(sg)
        else:
            self._is_wyckoff_list_fixed = False
            for sg in space_group_range:
                space_group = Sgroup(sg)
                if sg.wyckoff_preparation(self._nmpc):
                    self._space_group_range.append(sg)

        if len(self._space_group_range) == 0:
            raise ValueError(
                "No available space group matches the requirement of nmpc=%i,"
                " is_chiral=%s, wyckoff_list=%s, space_groups_allowed"
                "=%s" % (self._nmpc, str(self.is_chiral),
                         str(self._wyckoff_list),
                         str(self._space_groups_allowed)))
        print_time_log("Space group range from input: " +
                " ".join(map(str, self._space_group_range)))

    def get_space_group_randomly(self):
        sg = self._space_group_range[
            int(random.uniform(0, len(self._space_group_range)))]

        self._space_group = Sgroup(sg)
        self._space_group.wyckoff_preparation(self._nmpc)
        if not self._is_wyckoff_list_fixed:
            self._wyckoff_list = \
                self._space_group.wyckoff_selection(self._nmpc)

        self._space_group_selected_counter[sg] += 1
        return self._space_group

    def get_wyckoff_list_randomly(self):
        if not self._is_wyckoff_list_fixed:
            self._wyckoff_list = \
                self._space_group.wyckoff_selection(self._nmpc)

        return self._wyckoff_list

    def increment_space_group_counter(self):
        self._space_group_user_counter[
            self._space_group.space_group_number] += 1

class Sgroup():        
    def __init__(self,space_group_number=0,alternative_setting=0):
        '''
        elif space_group_number==:
            self.name=
            self.pgsym=['+x+y+z'.'-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0]]
            self.nwyc=
            self.wycs=['+x+y+z']
            self.wtran=[[0,0,0]]
            self.wmult=[]
            self.btype='P'
            self.blt=4
            self.nsg=
            self.sname=[]
            self.swycn=[]
            self.swyc_mult=[]
            self.snrep=
        elif space_group_number==:
            l=setting[alternative_setting]
            
            if l<0.5:         
                self.name=
                self.pgsym=['+x+y+z'.'-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0]]
                self.nwyc=
                self.wycs=['+x+y+z']
                self.wtran=[[0,0,0]]
                self.wmult=[]
            else:
                self.name=
                self.pgsym=['+x+y+z'.'-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0]]
                self.nwyc=
                self.wycs=['+x+y+z']
                self.wtran=[[0,0,0]]
                self.wmult=[]
            self.btype=''
            self.blt=4
            self.nsg=
            self.sname=[]
            self.swycn=[]
            self.swyc_mult=[]
            self.snrep=
        '''
        # TODO: Add the information for the rest of the space groups
        # TODO: Add alternative settings for space groups
        if space_group_number == -1:
            space_group_number = int(random.uniform(0,MAX_AVAILABLE_SPACE_GROUP+0.99999))
        
        self.space_group_number = space_group_number
        if space_group_number==0:
            self.name='random'
            self.pgsym=['+x+y+z']   #Rotation in symmetry operation
            self.trans=[[0,0,0]]    #Translation in symmetry operation
            self.nwyc=1             #No. of Wyckoff Positions
            self.wycs=['+x+y+z']    #Constraints on Wyckoff positions
            self.wtran=[[0,0,0]]    #Translation at Wyckoff positions
            self.wmult=[1]          #Multiplicity of Wyckoff positions
            self.nsg=1              #No. of Site symmetry group
            self.sname=['1']        #List of site symmetry name
            self.swycn=[1]          #Number of Wyckoff positions under the site symmetry group
            self.snrep=1            #The beginning when members in the site symmetry cannot be repeated
            self.swyc_mult=[1]      #Multiplicity of the Wyckoff positions under the space group (makes self.wmult obsolete)
            #If self.snrep>self.nsg-1, then all the Wyckoff positions under the site symmetry group can be repeated
            self.btype='P'          #Bravais centering type
            self.blt=-1             #Bravais Lattice type =-1 means random
        elif space_group_number==1:
            self.name='P1'
            self.pgsym=['+x+y+z']
            self.trans=[[0,0,0]]
            self.nwyc=1
            self.wycs=['+x+y+z']
            self.wtran=[[0,0,0]]
            self.wmult=[1]
            self.nsg=1              #No. of Site symmetry group
            self.sname=['1']        #List of site symmetry name
            self.swycn=[1]          #Number of Wyckoff positions under the site symmetry group
            self.swyc_mult=[1]
            self.snrep=1
            self.btype='P'
            self.blt=1
        elif space_group_number==2:
            self.name='P-1'
            self.pgsym=['+x+y+z','-x-y-z']
            self.trans=[[0,0,0],[0,0,0]]
            self.nwyc=9
            self.wycs=['+x+y+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[2,1,1,1,1,1,1,1,1]
            self.nsg=2
            self.sname=['1','-1']        
            self.swycn=[1,8]
            self.swyc_mult=[2,1]
            self.snrep=1
            self.btype='P'
            self.blt=1
        elif space_group_number==3:
#            l=setting[alternative_setting]
            l=0.25 #Only beta angle non-90
            if l<0.5:
                self.name='P2/b'
                self.pgsym=['+x+y+z','-x+y-z']
                self.trans=[[0,0,0],[0,0,0]]
                self.nwyc=5
                self.wycs=['+x+y+z','+0+y+0','+0+y+0','+0+y+0','+0+y+0']
                self.wtran=[[0,0,0],[0.5,0,0.5],[0.5,0,0],[0,0,0.5],[0,0,0]]
                self.wmult=[2,1,1,1,1]
                self.btype='P'
                self.blt=2
            else:
                self.name='P2/c'
                self.pgsym=['+x+y+z','-x-y+z']
                self.trans=[[0,0,0],[0,0,0]]
                self.nwyc=5
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z','+0+0+z']
                self.wtran=[[0,0,0],[0.5,0.5,0],[0,0.5,0],[0.5,0,0],[0,0,0]]
                self.wmult=[2,1,1,1,1]
                self.btype='P'
                self.blt=2
            self.nsg=2
            self.sname=['1','2']
            self.swycn=[1,4]
            self.swyc_mult=[2,1]
            self.snrep=1
        elif space_group_number==4:
#            l=setting[alternative_setting]
            l=0.25
            if l<0.5:                
                self.name='P2_1/b'
                self.pgsym=['+x+y+z','-x+y-z']
                self.trans=[[0,0,0],[0,0.5,0]]
                self.nwyc=1
                self.wycs=['+x+y+z']
                self.wtran=[[0,0,0]]
                self.wmult=[2]
                self.btype='P'
                self.blt=2
            if l>0.5:
                self.name='P2_1/c'
                self.pgsym=['+x+y+z','-x-y+z']
                self.trans=[[0,0,0],[0,0,0.5]]
                self.nwyc=1
                self.wycs=['+x+y+z']
                self.wtran=[[0,0,0]]
                self.wmult=[2]
                self.btype='P'
                self.blt=2 
            self.nsg=1
            self.sname=['1']
            self.swycn=[1]
            self.swyc_mult=[2]
            self.snrep=1
        elif space_group_number==5:
#            l=setting[alternative_setting]
            l=0.25
            if l<0.5:                
                self.name='C2/b'
                self.pgsym=['+x+y+z','-x+y-z']
                self.trans=[[0,0,0],[0,0,0]]
                self.nwyc=3
                self.wycs=['+x+y+z','+0+y+0','+0+y+0']
                self.wtran=[[0,0,0],[0,0,0.5],[0,0,0]]
                self.wmult=[2,1,1]
                self.btype='C'
                self.blt=2
            if l>0.5:
                self.name='P2_1/c'
                self.pgsym=['+x+y+z','-x+y-z']
                self.trans=[[0,0,0],[0,0,0]]
                self.nwyc=3
                self.wycs=['+x+y+z','+0+y+0','+0+y+0']
                self.wtran=[[0,0,0],[0,0,0.5],[0,0,0]]
                self.wmult=[2,1,1]
                self.btype='A'
                self.blt=2
            self.nsg=2
            self.sname=['1','2']
            self.swycn=[1,2]
            self.swyc_mult=[2,1]
            self.snrep=2
        elif space_group_number==6:
            l=0.25
            if l<0.5:                
                self.name='Pm/b'
                self.pgsym=['+x+y+z','+x-y+z']
                self.trans=[[0,0,0],[0,0,0]]
                self.nwyc=3
                self.wycs=['+x+y+z','+x+0+z','+x+0+z']
                self.wtran=[[0,0,0],[0,0.5,0],[0,0,0]]
                self.wmult=[2,1,1]
                self.btype='P'
                self.blt=2
            if l>0.5:
                self.name='Pm/c'
                self.pgsym=['+x+y+z','+x+y-z']
                self.trans=[[0,0,0],[0,0,0]]
                self.nwyc=3
                self.wycs=['+x+y+z','+x+y+0','+x+y+0']
                self.wtran=[[0,0,0],[0,0,0.5],[0,0,0]]
                self.wmult=[2,1,1]
                self.btype='P'
                self.blt=2   
            self.nsg=2
            self.sname=['1','m']
            self.swycn=[1,2]
            self.swyc_mult=[2,1]
            self.snrep=2
        elif space_group_number==7:
            l=0.25
            if l<0.5:                
                self.name='Pc/b'
                self.pgsym=['+x+y+z','+x-y+z']
                self.trans=[[0,0,0],[0,0,0.5]]
                self.nwyc=1
                self.wycs=['+x+y+z']
                self.wtran=[[0,0,0]]
                self.wmult=[2]
                self.btype='P'
                self.blt=2
            if l>0.5:
                self.name='Pc/c'
                self.pgsym=['+x+y+z','+x+y-z']
                self.trans=[[0,0,0],[0.5,0,0]]
                self.nwyc=1
                self.wycs=['+x+y+z']
                self.wtran=[[0,0,0]]
                self.wmult=[2]
                self.btype='P'
                self.blt=2
            self.nsg=1
            self.sname=['1']
            self.swycn=[1]
            self.swyc_mult=[2]
            self.snrep=1
        elif space_group_number==8:
            l=0.25
            if l<0.5:                
                self.name='Cm/b'
                self.pgsym=['+x+y+z','+x-y+z']
                self.trans=[[0,0,0],[0,0,0]]
                self.nwyc=2
                self.wycs=['+x+y+z','+x+0+z']
                self.wtran=[[0,0,0],[0,0,0]]
                self.wmult=[2,1]
                self.btype='C'
                self.blt=2
            if l>0.5:
                self.name='Cm/c'
                self.pgsym=['+x+y+z','+x+y-z']
                self.trans=[[0,0,0],[0,0,0]]
                self.nwyc=2
                self.wycs=['+x+y+z','+x+y+0']
                self.wtran=[[0,0,0],[0,0,0]]
                self.wmult=[2,1]
                self.btype='A'
                self.blt=2
            self.nsg=2
            self.sname=['1','m']
            self.swycn=[1,1]
            self.swyc_mult=[2,1]
            self.snrep=2
        elif space_group_number==9:
            l=0.25
            if l<0.5:                
                self.name='Cc/b'
                self.pgsym=['+x+y+z','+x-y+z']
                self.trans=[[0,0,0],[0,0,0.5]]
                self.nwyc=1
                self.wycs=['+x+y+z']
                self.wtran=[[0,0,0]]
                self.wmult=[2]
                self.btype='C'
                self.blt=2
            if l>0.5:
                self.name='Cc/c'
                self.pgsym=['+x+y+z','+x+y-z']
                self.trans=[[0,0,0],[0.5,0,0]]
                self.nwyc=1
                self.wycs=['+x+y+z']
                self.wtran=[[0,0,0]]
                self.wmult=[2]
                self.btype='A'
                self.blt=2
            self.nsg=1
            self.sname=['1']
            self.swycn=[1]
            self.swyc_mult=[2]
            self.snrep=1
        elif space_group_number==10:
            l=0.25
            if l<0.5:                
                self.name='P2/m/b'
                self.pgsym=['+x+y+z','-x+y-z','-x-y-z','+x-y+z']
                self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
                self.nwyc=15
                self.wycs=['+x+y+z','+x+0+z','+x+0+z','+0+y+0','+0+y+0','+0+y+0','+0+y+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0.5,0,0.5],[0,0,0.5],[0.5,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0],[0.5,0,0],[0,0,0.5],[0,0.5,0],[0,0,0]]
                self.wmult=[4,2,2,2,2,2,2,1,1,1,1,1,1,1,1]
                self.btype='P'
                self.blt=2
            if l>0.5:
                self.name='P2/m/c'
                self.pgsym=['+x+y+z','-x-y+z','-x-y-z','+x+y-z']
                self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
                self.nwyc=15
                self.wycs=['+x+y+z','+x+y+0','+x+y+0','+0+0+z','+0+0+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0.5],[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0],[0.5,0.5,0.5],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0],[0.5,0,0],[0,0,0.5],[0,0.5,0],[0,0,0]]
                self.wmult=[4,2,2,2,2,2,2,1,1,1,1,1,1,1,1]
                self.btype='P'
                self.blt=2    
            self.nsg=4
            self.sname=['1','m','2','2/m']
            self.swycn=[1,2,4,8]
            self.swyc_mult=[4,2,2,1]
            self.snrep=2
        elif space_group_number==11:
            l=0.25
            if l<0.5:                
                self.name='P2_1/m/b'
                self.pgsym=['+x+y+z','-x+y-z','-x-y-z','+x-y+z']
                self.trans=[[0,0,0],[0,0.5,0],[0,0,0],[0,0.5,0]]
                self.nwyc=6
                self.wycs=['+x+y+z','+x+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.25,0],[0.5,0,0.5],[0,0,0.5],[0.5,0,0],[0,0,0]]
                self.wmult=[4,2,2,2,2,2]
                self.btype='P'
                self.blt=2
            if l>0.5:
                self.name='P2_1/m/c'
                self.pgsym=['+x+y+z','-x-y+z','-x-y-z','+x+y-z']
                self.trans=[[0,0,0],[0,0,0.5],[0,0,0],[0,0,0.5]]
                self.nwyc=6
                self.wycs=['+x+y+z','+x+y+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0.25],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0]]
                self.wmult=[4,2,2,2,2,2]
                self.btype='P'
                self.blt=2    
            self.nsg=3
            self.sname=['1','m','-1']
            self.swycn=[1,1,4]
            self.swyc_mult=[4,2,2]
            self.snrep=2
        elif space_group_number==12:
            l=0.25
            if l<0.5:                
                self.name='C2/m/b'
                self.pgsym=['+x+y+z','-x+y-z','-x-y-z','+x-y+z']
                self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
                self.nwyc=10
                self.wycs=['+x+y+z','+x+0+z','+0+y+0','+0+y+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0],[0.25,0.25,0.5],[0.25,0.25,0],[0,0.5,0.5],[0,0,0.5],[0,0.5,0],[0,0,0]]
                self.wmult=[4,2,2,2,2,2,1,1,1,1]
                self.btype='C'
                self.blt=2
            if l>0.5:
                self.name='C2/m/c'
                self.pgsym=['+x+y+z','-x-y+z','-x-y-z','+x+y-z']
                self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
                self.nwyc=10
                self.wycs=['+x+y+z','+x+y+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0.5,0,0],[0,0,0],[0.5,0.25,0.25],[0,0.25,0.25],[0.5,0,0.5],[0.5,0,0],[0,0,0.5],[0,0,0]]
                self.wmult=[4,2,2,2,2,2,1,1,1,1]
                self.btype='A'
                self.blt=2    
            self.nsg=5
            self.sname=['1','m','2','-1','2/m']
            self.swycn=[1,1,2,2,4]
            self.swyc_mult=[4,2,2,2,1]
            self.snrep=3
        elif space_group_number==13:
            l=0.25
            if l<0.5:                
                self.name='P2/c/b'
                self.pgsym=['+x+y+z','-x+y-z','-x-y-z','+x-y+z']
                self.trans=[[0,0,0],[0,0,0.5],[0,0,0],[0,0,0.5]]
                self.nwyc=7
                self.wycs=['+x+y+z','+0+y+0','+0+y+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.5,0,0.25],[0,0,0.25],[0.5,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
                self.wmult=[4,2,2,2,2,2,2]
                self.btype='P'
                self.blt=2
            if l>0.5:
                self.name='P2/c/c'
                self.pgsym=['+x+y+z','-x-y+z','-x-y-z','+x+y-z']
                self.trans=[[0,0,0],[0.5,0,0],[0,0,0],[0.5,0,0]]
                self.nwyc=7
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.25,0.5,0],[0.25,0,0],[0,0.5,0],[0,0,0.5],[0,0.5,0.5],[0,0,0]]
                self.wmult=[4,2,2,2,2,2,2]
                self.btype='P'
                self.blt=2    
            self.nsg=3
            self.sname=['1','2','-1']
            self.swycn=[1,2,4]
            self.swyc_mult=[4,2,2]
            self.snrep=2
        elif space_group_number==14:
            l=0.25
            if l<0.5:                
                self.name='P2_1/c/b'
                self.pgsym=['+x+y+z','-x+y-z','-x-y-z','+x-y+z']
                self.trans=[[0,0,0],[0,0.5,0.5],[0,0,0],[0,0.5,0.5]]
                self.nwyc=5
                self.wycs=['+x+y+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.5,0,0.5],[0,0,0.5],[0.5,0,0],[0,0,0]]
                self.wmult=[4,2,2,2,2]
                self.btype='P'
                self.blt=2
            if l>=0.5: #Turn this back
                self.name='P2_1/c/c'
                self.pgsym=['+x+y+z','-x-y+z','-x-y-z','+x+y-z']
                self.trans=[[0,0,0],[0.5,0,0.5],[0,0,0],[0.5,0,0.5]]
                self.nwyc=5
                self.wycs=['+x+y+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0]]
                self.wmult=[4,2,2,2,2]
                self.btype='P'
                self.blt=2    
            self.nsg=2
            self.sname=['1','-1']
            self.swycn=[1,4]
            self.swyc_mult=[4,2]
            self.snrep=1
        elif space_group_number==15:
            l=0.25
            if l<0.5:                
                self.name='C2/c/b'
                self.pgsym=['+x+y+z','-x+y-z','-x-y-z','+x-y+z']
                self.trans=[[0,0,0],[0,0,0.5],[0,0,0],[0,0,0.5]]
                self.nwyc=6
                self.wycs=['+x+y+z','+0+y+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0.25],[0.25,0.25,0.5],[0.25,0.25,0],[0,0.5,0],[0,0,0]]
                self.wmult=[4,2,2,2,2,2]
                self.btype='C'
                self.blt=2
            if l>0.5:
                self.name='C2/c/c'
                self.pgsym=['+x+y+z','-x-y+z','-x-y-z','+x+y-z']
                self.trans=[[0,0,0],[0.5,0,0],[0,0,0],[0.5,0,0]]
                self.nwyc=6
                self.wycs=['+x+y+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.25,0,0],[0.5,0.25,0.25],[0,0.25,0.25],[0,0,0.5],[0,0,0]]
                self.wmult=[4,2,2,2,2,2]
                self.btype='A'
                self.blt=2    
            self.nsg=3
            self.sname=['1','2','-1']
            self.swycn=[1,1,4]
            self.swyc_mult=[4,2,2]
            self.snrep=2
        elif space_group_number==16:
            self.name='P222'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=21
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z','+0+0+z','+0+y+0','+0+y+0','+0+y+0','+0+y+0','+x+0+0','+x+0+0','+x+0+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0.5,0.5,0],[0,0.5,0],[0.5,0,0],[0,0,0],[0.5,0,0.5],[0.5,0,0],[0,0,0.5],[0,0,0],[0,0.5,0.5],[0,0.5,0],[0,0,0.5],[0,0,0],[0.5,0.5,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0,0,0.5],[0,0.5,0],[0.5,0,0],[0,0,0]]
            self.wmult=[4,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1]
            self.btype='P'
            self.blt=3    
            self.nsg=5
            self.sname=['1','..2','.2.','2..','222']
            self.swycn=[1,4,4,4,8]
            self.swyc_mult=[4,2,2,2,1]
            self.snrep=4
        elif space_group_number==17:
            self.name='P222_1'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0]]
            self.nwyc=5
            self.wycs=['+x+y+z','+0+y+0','+0+y+0','+x+0+0','+x+0+0']
            self.wtran=[[0,0,0],[0.5,0,0.25],[0,0,0.25],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2,2,2]
            self.btype='P'
            self.blt=3    
            self.nsg=3
            self.sname=['1','.2.','2..']
            self.swycn=[1,2,2]
            self.swyc_mult=[4,2,2]
            self.snrep=3
        elif space_group_number==18:
            self.name='P2_12_12'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2]
            self.btype='P'
            self.blt=3    
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,2]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==19:
            self.name='P2_12_12_1'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z']
            self.trans=[[0,0,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0]]
            self.nwyc=1
            self.wycs=['+x+y+z']
            self.wtran=[[0,0,0]]
            self.wmult=[4]
            self.btype='P'
            self.blt=3    
            self.nsg=1
            self.sname=['1']
            self.swycn=[1]
            self.swyc_mult=[4]
            self.snrep=1
        elif space_group_number==20:
            self.name='C222_1'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+y+0','+x+0+0']
            self.wtran=[[0,0,0],[0,0,0.25],[0,0,0]]
            self.wmult=[4,2,2]
            self.btype='C'
            self.blt=3    
            self.nsg=3
            self.sname=['1','.2.','2..']
            self.swycn=[1,1,1]
            self.swyc_mult=[4,2,2]
            self.snrep=3
        elif space_group_number==21:
            self.name='C222'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=12
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z','+0+y+0','+0+y+0','+x+0+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0.25,0.25,0],[0,0.5,0],[0,0,0],[0,0,0.5],[0,0,0],[0,0,0.5],[0,0,0],[0,0,0.5],[0.5,0,0.5],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2,2,2,2,2,2,1,1,1,1]
            self.btype='C'
            self.blt=3    
            self.nsg=5
            self.sname=['1','..2','.2.','2..','222']
            self.swycn=[1,3,2,2,4]
            self.swyc_mult=[4,2,2,2,1]
            self.snrep=4
        elif space_group_number==22:
            self.name='F222'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=11
            self.wycs=['+x+y+z','+x+0+0','+x+0+0','+0+y+0','+0+y+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.25,0.25],[0,0,0],[0.25,0,0.25],[0,0,0],[0.25,0.25,0],[0,0,0],[0.25,0.25,0.75],[0.25,0.25,0.25],[0,0,0.5],[0,0,0]]
            self.wmult=[4,2,2,2,2,2,2,1,1,1,1]
            self.btype='F'
            self.blt=3    
            self.nsg=5
            self.sname=['1','2..','.2.','..2','222']
            self.swycn=[1,2,2,2,4]
            self.swyc_mult=[4,2,2,2,1]
            self.snrep=4
        elif space_group_number==23:
            self.name='I222'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=11
            self.wycs=['+x+y+z','+x+0+0','+x+0+0','+0+y+0','+0+y+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.5],[0,0,0],[0.5,0,0],[0,0,0],[0,0.5,0],[0,0,0],[0,0.5,0],[0,0,0.5],[0.5,0,0],[0,0,0]]
            self.wmult=[4,2,2,2,2,2,2,1,1,1,1]
            self.btype='I'
            self.blt=3    
            self.nsg=5
            self.sname=['1','2..','.2.','..2','222']
            self.swycn=[1,2,2,2,4]
            self.swyc_mult=[4,2,2,2,1]
            self.snrep=4
        elif space_group_number==24:
            self.name='I2_12_12_1'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z']
            self.trans=[[0,0,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0]]
            self.nwyc=4
            self.wycs=['+x+y+z','+0+0+z','+0+y+0','+x+0+0']
            self.wtran=[[0,0,0],[0,0.25,0],[0.25,0,0],[0,0,0.25]]
            self.wmult=[4,2,2,2]
            self.btype='I'
            self.blt=3    
            self.nsg=4
            self.sname=['1','..2','.2.','2..']
            self.swycn=[1,1,1,1]
            self.swyc_mult=[4,2,2,2]
            self.snrep=4
        elif space_group_number==25:
            self.name='Pmm2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=9
            self.wycs=['+x+y+z','+0+y+z','+0+y+z','+x+0+z','+x+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0.5,0,0],[0,0,0],[0,0.5,0],[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2,2,2,1,1,1,1]
            self.btype='P'
            self.blt=3    
            self.nsg=4
            self.sname=['1','m..','.m.','mm2']
            self.swycn=[1,2,2,4]
            self.swyc_mult=[4,2,2,1]
            self.snrep=4
        elif space_group_number==26: #start test here
            self.name='Pmc2_1'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+y+z','+0+y+z']
            self.wtran=[[0,0,0],[0.5,0,0],[0,0,0]]
            self.wmult=[4,2,2]
            self.btype='P'
            self.blt=3    
            self.nsg=2
            self.sname=['1','m..']
            self.swycn=[1,2]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==27:
            self.name='Pcc2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5]]
            self.nwyc=5
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2,2,2]
            self.btype='P'
            self.blt=3    
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,4]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==28:
            self.name='Pma2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0,0],[0.5,0,0]]
            self.nwyc=4
            self.wycs=['+x+y+z','+0+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0.25,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2,2]
            self.btype='P'
            self.blt=3    
            self.nsg=3
            self.sname=['1','m..','..2']
            self.swycn=[1,1,2]
            self.swyc_mult=[4,2,2]
            self.snrep=3 
        elif space_group_number==29:
            self.name='Pca2_1'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0.5],[0.5,0,0],[0.5,0,0.5]]
            self.nwyc=1
            self.wycs=['+x+y+z']
            self.wtran=[[0,0,0]]
            self.wmult=[4]
            self.btype='P'
            self.blt=3    
            self.nsg=1
            self.sname=['1']
            self.swycn=[1]
            self.swyc_mult=[4]
            self.snrep=1 
        elif space_group_number==30:
            self.name='Pnc2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0.5,0.5],[0,0.5,0.5]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0.5,0,0],[0,0,0]]
            self.wmult=[4,2,2]
            self.btype='P'
            self.blt=3    
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,2]
            self.swyc_mult=[4,2]
            self.snrep=3
        elif space_group_number==31:
            self.name='Pmn2_1'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0,0.5],[0.5,0,0.5],[0,0,0]]
            self.nwyc=2
            self.wycs=['+x+y+z','+0+y+z']
            self.wtran=[[0,0,0],[0,0,0]]
            self.wmult=[4,2]
            self.btype='P'
            self.blt=3    
            self.nsg=2
            self.sname=['1','m..']
            self.swycn=[1,1]
            self.swyc_mult=[4,2]
            self.snrep=3
        elif space_group_number==32:
            self.name='Pba2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2]
            self.btype='P'
            self.blt=3    
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,2]
            self.swyc_mult=[4,2]
            self.snrep=3
        elif space_group_number==33:
            self.name='Pna2_1'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0.5],[0.5,0.5,0],[0.5,0.5,0.5]]
            self.nwyc=1
            self.wycs=['+x+y+z']
            self.wtran=[[0,0,0]]
            self.wmult=[4]
            self.btype='P'
            self.blt=3    
            self.nsg=1
            self.sname=['1']
            self.swycn=[1]
            self.swyc_mult=[4]
            self.snrep=1
        elif space_group_number==34:
            self.name='Pnn2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2]
            self.btype='P'
            self.blt=3    
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,2]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==35:
            self.name='Cmm2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=6
            self.wycs=['+x+y+z','+0+y+z','+x+0+z','+0+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0.25,0.25,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2,2,1,1]
            self.btype='C'
            self.blt=3    
            self.nsg=5
            self.sname=['1','m..','.m.','..2','mm2']
            self.swycn=[1,1,1,1,2]
            self.swyc_mult=[4,2,2,2,1]
            self.snrep=5       
        elif space_group_number==36:
            self.name='Cmc2_1'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0]]
            self.nwyc=2
            self.wycs=['+x+y+z','+0+y+z']
            self.wtran=[[0,0,0],[0,0,0]]
            self.wmult=[4,2]
            self.btype='C'
            self.blt=3    
            self.nsg=2
            self.sname=['1','m..']
            self.swycn=[1,1]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==37:
            self.name='Ccc2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5]]
            self.nwyc=4
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0.25,0.25,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2,2]
            self.btype='C'
            self.blt=3    
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,3]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==38:
            self.name='Amm2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=6
            self.wycs=['+x+y+z','+0+y+z','+0+y+z','+x+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0.5,0,0],[0,0,0],[0,0,0],[0.5,0,0],[0,0,0]]
            self.wmult=[4,2,2,2,1,1]
            self.btype='A'
            self.blt=3    
            self.nsg=4
            self.sname=['1','m..','.m.','mm2']
            self.swycn=[1,2,1,2]
            self.swyc_mult=[4,2,2,1]
            self.snrep=4
        elif space_group_number==39:
            self.name='Aem2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0.5,0],[0,0.5,0]]
            self.nwyc=4
            self.wycs=['+x+y+z','+x+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.25,0],[0.5,0,0],[0,0,0]]
            self.wmult=[4,2,2,2]
            self.btype='A'
            self.blt=3    
            self.nsg=3
            self.sname=['1','.m.','..2']
            self.swycn=[1,1,2]
            self.swyc_mult=[4,2,2]
            self.snrep=3
        elif space_group_number==40:
            self.name='Ama2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0,0],[0.5,0,0]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+y+z','+0+0+z']
            self.wtran=[[0,0,0],[0.25,0,0],[0,0,0]]
            self.wmult=[4,2,2]
            self.btype='A'
            self.blt=3    
            self.nsg=3
            self.sname=['1','m..','..2']
            self.swycn=[1,1,1]
            self.swyc_mult=[4,2,2]
            self.snrep=3
        elif space_group_number==41:
            self.name='Aea2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=2
            self.wycs=['+x+y+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0]]
            self.wmult=[4,2]
            self.btype='A'
            self.blt=3    
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,1]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==42:
            self.name='Fmm2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=5
            self.wycs=['+x+y+z','+x+0+z','+0+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0.25,0.25,0],[0,0,0]]
            self.wmult=[4,2,2,2,1]
            self.btype='F'
            self.blt=3    
            self.nsg=5
            self.sname=['1','.m.','m..','..2','mm2']
            self.swycn=[1,1,1,1,1]
            self.swyc_mult=[4,2,2,2,1]
            self.snrep=5
        elif space_group_number==43:
            self.name='Fdd2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.25,0.25,0.25],[0.25,0.25,0.25]]
            self.nwyc=2
            self.wycs=['+x+y+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0]]
            self.wmult=[4,2]
            self.btype='F'
            self.blt=3    
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,1]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==44:
            self.name='Imm2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=5
            self.wycs=['+x+y+z','+0+y+z','+x+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2,1,1]
            self.btype='I'
            self.blt=3    
            self.nsg=4
            self.sname=['1','m..','.m.','mm2']
            self.swycn=[1,1,1,2]
            self.swyc_mult=[4,2,2,1]
            self.snrep=4
        elif space_group_number==45:
            self.name='Iba2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,2]
            self.btype='I'
            self.blt=3    
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,2]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==46:
            self.name='Ima2'
            self.pgsym=['+x+y+z','-x-y+z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0,0],[0.5,0,0]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+y+z','+0+0+z']
            self.wtran=[[0,0,0],[0.25,0,0],[0,0,0]]
            self.wmult=[4,2,2]
            self.btype='I'
            self.blt=3    
            self.nsg=3
            self.sname=['1','m..','..2']
            self.swycn=[1,1,1]
            self.swyc_mult=[4,2,2]
            self.snrep=3
        elif space_group_number==47:
            self.name='Pmmm'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=27
            self.wycs=['+x+y+z','+x+y+0','+x+y+0','+x+0+z','+x+0+z','+0+y+z','+0+y+z',
                       '+0+0+z','+0+0+z','+0+0+z','+0+0+z',
                       '+0+y+0','+0+y+0','+0+y+0','+0+y+0',
                       '+x+0+0','+x+0+0','+x+0+0','+x+0+0',
                       '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.5],[0,0,0],[0,0.5,0],[0,0,0],[0.5,0,0],[0,0,0],
                        [0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0],
                        [0.5,0,0.5],[0.5,0,0],[0,0,0.5],[0,0,0],
                        [0,0.5,0.5],[0,0.5,0],[0,0,0.5],[0,0,0],
                        [0.5,0.5,0.5],[0,0.5,0.5],[0.5,0.5,0],[0,0.5,0],
                        [0.5,0,0.5],[0,0,0.5],[0.5,0,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,4,4,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1]
            self.btype='P'
            self.blt=3
            self.nsg=8
            self.sname=['1','..m','.m.','m..','mm2','m2m','2mm','mmm']
            self.swycn=[1,2,2,2,4,4,4,8]
            self.swyc_mult=[8,4,4,4,2,2,2,1]
            self.snrep=7
        elif space_group_number==48:
            l=setting[alternative_setting]
            if l<0.5:                
                self.name='Pnnn_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]
                self.nwyc=13
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+y+0','+0+y+0','+x+0+0','+x+0+0',
                           '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0.5,0,0],[0,0,0],[0,0,0.5],[0,0,0],
                            [0.75,0.75,0.75],[0.25,0.25,0.25],[0,0.5,0],[0,0,0.5],[0.5,0,0],[0,0,0]]
                self.wmult=[8,4,4,4,4,4,4,4,4,2,2,2,2]
                self.btype='P'
                self.blt=3
            if l>0.5:
                self.name='Pnnn_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]]
                self.nwyc=13
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+y+0','+0+y+0','+x+0+0','+x+0+0',
                           '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.25,0.75,0],[0.25,0.25,0],[0.75,0,0.25],[0.25,0,0.25],[0,0.25,0.75],[0,0.25,0.25],
                            [0,0,0],[0.5,0.5,0.5],[0.25,0.75,0.25],[0.25,0.25,0.75],[0.75,0.25,0.25],[0.25,0.25,0.25]]
                self.wmult=[8,4,4,4,4,4,4,4,4,2,2,2,2]
                self.btype='P'
                self.blt=3    
            self.nsg=6
            self.sname=['1','..2','.2.','2..','-1','222']
            self.swycn=[1,2,2,2,2,4]
            self.swyc_mult=[8,4,4,4,4,2]
            self.snrep=4
        elif space_group_number==49:
            self.name='Pccm'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5]]
            self.nwyc=18
            self.wycs=['+x+y+z','+x+y+0',
                       '+0+0+z','+0+0+z','+0+0+z','+0+0+z',
                       '+0+y+0','+0+y+0',
                       '+x+0+0','+x+0+0',
                       '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0.5,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0],
                        [0.5,0,0.25],[0,0,0.25],[0,0.5,0.25],[0,0,0.25],
                        [0.5,0.5,0.25],[0,0.5,0.25],[0.5,0,0.25],[0,0,0.25],
                        [0.5,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2]
            self.btype='P'
            self.blt=3
            self.nsg=7
            self.sname=['1','..m','..2','.2.','2..','222','..2/m']
            self.swycn=[1,1,4,2,2,4,4]
            self.swyc_mult=[8,4,4,4,4,2,2]
            self.snrep=5
        elif space_group_number==50:
            l=setting[alternative_setting]
            if l<0.5:                
                self.name='Pban_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0]]
                self.nwyc=13
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+y+0','+0+y+0','+x+0+0','+x+0+0',
                           '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0,0.5],[0,0,0],[0,0,0.5],[0,0,0],
                            [0.25,0.25,0.5],[0.25,0.25,0],[0,0,0.5],[0.5,0,0.5],[0.5,0,0],[0,0,0]]
                self.wmult=[8,4,4,4,4,4,4,4,4,2,2,2,2]
                self.btype='P'
                self.blt=3
            if l>0.5:
                self.name='Pban_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0]]
                self.nwyc=13
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+y+0','+0+y+0','+x+0+0','+x+0+0',
                           '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.25,0.75,0],[0.25,0.25,0],[0.25,0,0.5],[0.25,0,0],[0,0.25,0.5],[0,0.25,0],
                            [0,0,0.5],[0,0,0],[0.25,0.25,0.5],[0.75,0.25,0.5],[0.75,0.25,0],[0.25,0.25,0]]
                self.wmult=[8,4,4,4,4,4,4,4,4,2,2,2,2]
                self.btype='P'
                self.blt=3    
            self.nsg=6
            self.sname=['1','..2','.2.','2..','-1','222']
            self.swycn=[1,2,2,2,2,4]
            self.swyc_mult=[8,4,4,4,4,2]
            self.snrep=4
        elif space_group_number==51:
            self.name='Pmma'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0,0],[0,0,0],[0.5,0,0],[0,0,0],[0.5,0,0],[0,0,0],[0.5,0,0]]
            self.nwyc=12
            self.wycs=['+x+y+z','+0+y+z','+x+0+z','+x+0+z','+0+y+0','+0+y+0',
                       '+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0.25,0,0],[0,0.5,0],[0,0,0],[0,0,0.5],[0,0,0],
                        [0.25,0.5,0],[0.25,0,0],[0,0.5,0.5],[0,0,0.5],[0,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,4,2,2,2,2,2,2]
            self.btype='P'
            self.blt=3
            self.nsg=6
            self.sname=['1','m..','.m.','.2.','mm2','.2/m.']
            self.swycn=[1,1,2,2,2,4]
            self.swyc_mult=[8,4,4,4,2,2]
            self.snrep=5
        elif space_group_number==52:
            self.name='Pnna'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0,0],[0.5,0.5,0.5],[0,0.5,0.5],[0,0,0],[0.5,0,0],[0.5,0.5,0.5],[0,0.5,0.5]]
            self.nwyc=5
            self.wycs=['+x+y+z','+x+0+0','+0+0+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.25,0.25],[0.25,0,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4]
            self.btype='P'
            self.blt=3
            self.nsg=4
            self.sname=['1','2..','..2','-1']
            self.swycn=[1,1,1,2]
            self.swyc_mult=[8,4,4,4]
            self.snrep=3
        elif space_group_number==53:
            self.name='Pmna'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0,0.5],[0.5,0,0.5],[0,0,0],[0,0,0],[0.5,0,0.5],[0.5,0,0.5],[0,0,0]]
            self.nwyc=9
            self.wycs=['+x+y+z','+0+y+z','+0+y+0','+x+0+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0.25,0,0.25],[0,0.5,0],[0,0,0],[0,0.5,0],[0,0.5,0.5],[0.5,0,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=3
            self.nsg=5
            self.sname=['1','m..','.2.','2..','2/m..']
            self.swycn=[1,1,1,2,4]
            self.swyc_mult=[8,4,4,4,2]
            self.snrep=4
        elif space_group_number==54:
            self.name='Pcca'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0,0],[0,0,0.5],[0.5,0,0.5],[0,0,0],[0.5,0,0],[0,0,0.5],[0.5,0,0.5]]
            self.nwyc=6
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+y+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0.25,0.5,0],[0.25,0,0],[0,0,0.25],[0,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,4]
            self.btype='P'
            self.blt=3
            self.nsg=4
            self.sname=['1','..2','.2.','-1']
            self.swycn=[1,2,1,2]
            self.swyc_mult=[8,4,4,4]
            self.snrep=3
        elif space_group_number==55:
            self.name='Pbam'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=9
            self.wycs=['+x+y+z','+x+y+0','+x+y+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.5],[0,0,0],[0,0.5,0],[0,0,0],[0,0.5,0.5],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=3
            self.nsg=4
            self.sname=['1','..m','..2','..2/m']
            self.swycn=[1,2,2,4]
            self.swyc_mult=[8,4,4,2]
            self.snrep=3
        elif space_group_number==56:
            self.name='Pccn'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5],[0,0,0],[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]]
            self.nwyc=5
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0.25,0.75,0],[0.25,0.25,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4]
            self.btype='P'
            self.blt=3
            self.nsg=3
            self.sname=['1','..2','-1']
            self.swycn=[1,2,2]
            self.swyc_mult=[8,4,4]
            self.snrep=2
        elif space_group_number==57:
            self.name='Pbcm'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0.5,0.5],[0,0.5,0],[0,0,0],[0,0,0.5],[0,0.5,0.5],[0,0.5,0]]
            self.nwyc=5
            self.wycs=['+x+y+z','+x+y+0','+x+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.25],[0,0.25,0],[0.5,0,0],[0,0,0]]
            self.wmult=[8,4,4,4,4]
            self.btype='P'
            self.blt=3
            self.nsg=4
            self.sname=['1','..m','2..','-1']
            self.swycn=[1,1,1,2]
            self.swyc_mult=[8,4,4,4]
            self.snrep=3
        elif space_group_number==58:
            self.name='Pnnm'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5]]
            self.nwyc=8
            self.wycs=['+x+y+z','+x+y+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0],[0,0,0],[0,0.5,0.5],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=3
            self.nsg=4
            self.sname=['1','..m','..2','..2/m']
            self.swycn=[1,1,2,4]
            self.swyc_mult=[8,4,4,2]
            self.snrep=3
        elif space_group_number==59:
            l=setting[alternative_setting]
            if l<0.5:                
                self.name='Pmmn_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0]]
                self.nwyc=7
                self.wycs=['+x+y+z','+x+0+z','+0+y+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0.25,0.25,0.5],[0.25,0.25,0],[0,0.5,0],[0,0,0]]
                self.wmult=[8,4,4,4,4,2,2]
                self.btype='P'
                self.blt=3
            if l>0.5:
                self.name='Pmmn_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0,0.5,0],[0.5,0,0],[0,0,0],[0.5,0.5,0],[0,0.5,0],[0.5,0,0]]
                self.nwyc=7
                self.wycs=['+x+y+z','+x+0+z','+0+y+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.25,0],[0.25,0,0],[0,0,0.5],[0,0,0],[0.25,0.75,0],[0.25,0.25,0]]
                self.wmult=[8,4,4,4,4,2,2]
                self.btype='P'
                self.blt=3    
            self.nsg=5
            self.sname=['1','.m.','m..','-1','mm2']
            self.swycn=[1,1,1,2,2]
            self.swyc_mult=[8,4,4,4,2]
            self.snrep=3
        elif space_group_number==60:
            self.name='Pbcn'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0.5,0.5],[0,0,0.5],[0.5,0.5,0],[0,0,0],[0.5,0.5,0.5],[0,0,0.5],[0.5,0.5,0]]
            self.nwyc=4
            self.wycs=['+x+y+z','+0+y+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.25],[0,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,4]
            self.btype='P'
            self.blt=3
            self.nsg=3
            self.sname=['1','.2.','-1']
            self.swycn=[1,1,2]
            self.swyc_mult=[8,4,4]
            self.snrep=2
        elif space_group_number==61:
            self.name='Pbca'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0],[0,0,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4]
            self.btype='P'
            self.blt=3
            self.nsg=2
            self.sname=['1','-1']
            self.swycn=[1,2]
            self.swyc_mult=[8,4]
            self.snrep=1
        elif space_group_number==62:
            self.name='Pnma'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0,0.5],[0,0.5,0],[0.5,0.5,0.5],[0,0,0],[0.5,0,0.5],[0,0.5,0],[0.5,0.5,0.5]]
            self.nwyc=4
            self.wycs=['+x+u+z','+x+0+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.25,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4]
            self.btype='P'
            self.blt=3
            self.nsg=3
            self.sname=['1','.m.','-1']
            self.swycn=[1,1,2]
            self.swyc_mult=[8,4,4]
            self.snrep=2
        elif space_group_number==63:
            self.name='Cmcm'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0]]
            self.nwyc=8
            self.wycs=['+x+y+z','+x+y+0','+0+y+z','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.25],[0,0,0],[0,0,0],[0.25,0.25,0],[0,0,0.25],[0,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2]
            self.btype='C'
            self.blt=3
            self.nsg=7
            self.sname=['1','..m','m..','2..','-1','m2m','2/m..']
            self.swycn=[1,1,1,1,1,1,2]
            self.swyc_mult=[8,4,4,4,4,2,2]
            self.snrep=6
        elif space_group_number==64:
            self.name='Cmce'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0.5,0.5],[0,0.5,0.5],[0,0,0],[0,0,0],[0,0.5,0.5],[0,0.5,0.5],[0,0,0]]
            self.nwyc=7
            self.wycs=['+x+y+z','+0+y+z','+0+y+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0.25,0,0.25],[0,0,0],[0.25,0.25,0],[0.5,0,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2]
            self.btype='C'
            self.blt=3
            self.nsg=6
            self.sname=['1','m..','.2.','2..','-1','2/m..']
            self.swycn=[1,1,1,1,1,2]
            self.swyc_mult=[8,4,4,4,4,2]
            self.snrep=4
        elif space_group_number==65:
            self.name='Cmmm'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=18
            self.wycs=['+x+y+z','+x+y+0','+x+y+0','+x+0+z','+0+y+z',
                       '+0+0+z','+0+0+z','+0+0+z',
                       '+0+y+0','+0+y+0','+x+0+0','+x+0+0',
                       '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0],
                        [0.25,0.25,0],[0,0.5,0],[0,0,0],
                        [0,0,0.5],[0,0,0],[0,0,0.5],[0,0,0],
                        [0.25,0.25,0.5],[0.25,0.25,0],[0,0,0.5],[0.5,0,0.5],[0.5,0,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,4,2,2,2,2,2,2,2,2,1,1,1,1]
            self.btype='C'
            self.blt=3
            self.nsg=10
            self.sname=['1','..m','.m.','m..','..2','mm2','m2m','2mm','..2/m','mmm']
            self.swycn=[1,2,1,1,1,2,2,2,2,4]
            self.swyc_mult=[8,4,4,4,4,2,2,2,2,1]
            self.snrep=8
        elif space_group_number==66:
            self.name='Cccm'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5]]
            self.nwyc=13
            self.wycs=['+x+y+z','+x+y+0','+0+0+z','+0+0+z','+0+0+z','+0+y+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0.25,0.25,0],[0,0.5,0],[0,0,0],[0,0,0.25],[0,0,0.25],
                        [0.25,0.75,0],[0.25,0.25,0],[0,0.5,0],[0,0,0],[0,0.5,0.25],[0,0,0.25]]
            self.wmult=[8,4,4,4,4,4,4,2,2,2,2,2,2]
            self.btype='C'
            self.blt=3
            self.nsg=7
            self.sname=['1','..m','..2','.2.','2..','..2/m','222']
            self.swycn=[1,1,3,1,1,4,2]
            self.swyc_mult=[8,4,4,4,4,2,2]
            self.snrep=5
        elif space_group_number==67:
            self.name='Cmm3'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0.5,0],[0,0.5,0],[0,0,0],[0,0,0],[0,0.5,0],[0,0.5,0],[0,0,0]]
            self.nwyc=15
            self.wycs=['+x+y+z','+x+0+z','+0+y+z','+0+0+z','+0+y+0','+0+y+0','+x+0+0','+x+0+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.25,0],[0,0,0],[0.25,0,0],[0.25,0,0.5],[0.25,0,0],[0,0,0.5],[0,0,0],[0,0.25,0],[0.25,0.25,0.5],[0.25,0.25,0],[0,0,0.5],[0,0,0],[0.25,0,0.5],[0.25,0,0]]
            self.wmult=[8,4,4,4,4,4,4,4,2,2,2,2,2,2,2]
            self.btype='C'
            self.blt=3
            self.nsg=10
            self.sname=['1','.m.','m..','..2','.2.','2..','mm2','.2/m.','2/m..','222']
            self.swycn=[1,1,1,1,2,2,1,2,2,2]
            self.swyc_mult=[8,4,4,4,4,4,2,2,2,2]
            self.snrep=7
        elif space_group_number==68:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='Ccce_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0,0,0],[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5],[0,0.5,0.5],[0.5,0,0.5]]
                self.nwyc=9
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+y+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.25,0.25,0],[0,0,0],[0,0,0],[0,0,0],[0,0.25,0.25],[0.25,0,0.25],[0,0,0.5],[0,0,0]]
                self.wmult=[8,4,4,4,4,4,4,2,2]
            else:
                self.name='Ccce_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0.5,0,0,],[0,0,0.5],[0.5,0,0.5],[0,0,0],[0.5,0,0],[0,0,0.5],[0.5,0,0.5]] 
                self.nwyc=9
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+y+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.25,0,0],[0,0.25,0],[0,0,0.25],[0,0.25,0.25],[0,0,0],[0.25,0.75,0],[0,0.25,0.75],[0,0.25,0.25]]
                self.wmult=[8,4,4,4,4,4,4,2,2]
            self.btype='C'
            self.blt=3
            self.nsg=6
            self.sname=['1','..2','.2.','2..','-1','222']
            self.swycn=[1,2,1,1,2,2]
            self.swyc_mult=[8,4,4,4,4,2]
            self.snrep=4
        elif space_group_number==69:
            self.name="Fmmm"
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=16
            self.wycs=['+x+y+z','+x+y+0','+x+0+z','+0+y+z','+x+0+0','+0+y+0','+0+0+z','+0+0+z','+0+y+0','+x+0+0',
                       '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0.25,0.25],[0.25,0,0.25,],[0.25,0.25,0],[0,0,0],[0,0,0],[0,0,0],
                        [0.25,0.25,0.25],[0.25,0.25,0],[0.25,0,0.25],[0,0.25,0.25],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,4,4,2,2,2,2,2,2,2,1,1]
            self.btype='F'
            self.blt=3
            self.nsg=15
            self.sname=['1','..m','.m.','m..','2..','.2.','..2','mm2','m2m','2mm','222','..2/m','.2/m','2/m..','mmm']
            self.swycn=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,2]
            self.swyc_mult=[8,4,4,4,4,4,4,2,2,2,2,2,2,2,1]
            self.snrep=10
        elif space_group_number==70:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='Fddd_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.25,0.25,0.25],[0.25,0.25,0.25],[0.25,0.25,0.25],[0.25,0.25,0.25]]
                self.nwyc=8
                self.wycs=['+x+y+z','+0+0+z','+0+y+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.625,0.625,0.625],[0.125,0.125,0.125],[0,0,0.5],[0,0,0]]
                self.wmult=[8,4,4,4,4,4,2,2]
            else:
                self.name='Fddd_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
                self.trans=[[0,0,0],[0.75,0.75,0],[0.75,0,0.75],[0,0.75,0.75],[0,0,0],[0.25,0.25,0],[0.25,0,0.25],[0,0.25,0.25]] 
                self.nwyc=8
                self.wycs=['+x+y+z','+0+0+z','+0+y+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.125,0.125,0],[0.125,0,0.125],[0,0.125,0.125],[0.5,0.5,0.5],[0,0,0],[0.125,0.125,0.625],[0.125,0.125,0.125]]
                self.wmult=[8,4,4,4,4,4,2,2]
            self.btype='F'
            self.blt=3
            self.nsg=6
            self.sname=['1','..2','.2.','2..','-1','222']
            self.swycn=[1,1,1,1,2,2]
            self.swyc_mult=[8,4,4,4,4,2]
            self.snrep=4
        elif space_group_number==71:
            self.name='Immm'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=15
            self.wycs=['+x+y+z','+x+y+0','+x+0+z','+0+y+z','+0+0+0',
                       '+0+0+z','+0+0+z','+0+y+0','+0+y+0','+x+0+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.25,0.25,0.25],[0.5,0,0],[0,0,0],[0,0,0.5],[0,0,0],[0,0.5,0],[0,0,0],[0.5,0,0.5],[0.5,0.5,0],[0,0.5,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2,2,2,2,1,1,1,1]
            self.btype='I'
            self.blt=3
            self.nsg=9
            self.sname=['1','..m','.m.','m..','-1','mm2','m2m','2mm','mmm']
            self.swycn=[1,1,1,1,1,2,2,2,4]
            self.swyc_mult=[8,4,4,4,4,2,2,2,1]
            self.snrep=8
        elif space_group_number==72:
            self.name='Ibam'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=11
            self.wycs=['+x+u+z','+x+y+0','+0+0+z','+0+0+z','+0+y+0','+x+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0],[0,0,0],[0,0,0.25],[0,0,0.25],[0.25,0.25,0.25],[0.5,0,0],[0,0,0],[0.5,0.,0.25],[0,0,0.25]]
            self.wmult=[8,4,4,4,4,4,4,2,2,2,2]
            self.btype='I'
            self.blt=3
            self.nsg=8
            self.sname=['1','..m','..2','.2.','2..','-1','..2/m','222']
            self.swycn=[1,1,2,1,1,1,2,2]
            self.swyc_mult=[8,4,4,4,4,4,2,2]
            self.snrep=5
        elif space_group_number==73:
            self.name='Ibca'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0],[0,0,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0]]
            self.nwyc=6
            self.wycs=['+x+u+z','+0+0+z','+0+y+0','+x+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.25,0],[0.25,0,0],[0,0,0.25],[0.25,0.25,0.25],[0,0,0]]
            self.wmult=[8,4,4,4,4,4]
            self.btype='I'
            self.blt=3
            self.nsg=5
            self.sname=['1','..2','.2.','2..','-1']
            self.swycn=[1,1,1,1,2]
            self.swyc_mult=[8,4,4,4,4]
            self.snrep=4
        elif space_group_number==74:
            self.name='Imma'
            self.pgsym=['+x+y+z','-x-y+z','-x+y-z','+x-y-z','-x-y-z','+x+y-z','+x-y+z','-x+y+z']
            self.trans=[[0,0,0],[0,0.5,0],[0,0.5,0],[0,0,0],[0,0,0],[0,0.5,0],[0,0.5,0],[0,0,0]]
            self.nwyc=10
            self.wycs=['+x+y+z','+x+0+z','+0+y+z','+0+y+0','+x+0+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.25,0],[0,0,0],[0.25,0,0.25],[0,0,0],[0,0.25,0],[0.25,0.25,0.75],[0.25,0.25,0.25],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2,2,2]
            self.btype='I'
            self.blt=3
            self.nsg=8
            self.sname=['1','.m.','m..','.2.','2..','mm2','.2/m.','2/m..']
            self.swycn=[1,1,1,1,1,1,2,2]
            self.swyc_mult=[8,4,4,4,4,2,2,2]
            self.snrep=6
        elif space_group_number==75:
            self.name='P4'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=4
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[4,2,1,1]
            self.btype='P'
            self.blt=4
            self.nsg=3
            self.sname=['1','2..','4..']
            self.swycn=[1,1,2]
            self.swyc_mult=[4,2,1]
            self.snrep=3
        elif space_group_number==76:
            self.name='P4_1'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0,0.25],[0,0,0.75]]
            self.nwyc=1
            self.wycs=['+x+y+z']
            self.wtran=[[0,0,0]]
            self.wmult=[4]
            self.btype='P'
            self.blt=4
            self.nsg=1
            self.sname=['1']
            self.swycn=[1]
            self.swyc_mult=[4]
            self.snrep=1
        elif space_group_number==77:
            self.name='P4_2'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5]]
            self.nwyc=4
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[4,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=2
            self.sname=['1','2..']
            self.swycn=[1,3]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==78:
            self.name='P4_3'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0,0.75],[0,0,0.25]]
            self.nwyc=1
            self.wycs=['+x+y+z']
            self.wtran=[[0,0,0]]
            self.wmult=[4]
            self.btype='P'
            self.blt=4
            self.nsg=1
            self.sname=['1']
            self.swycn=[1]
            self.swyc_mult=[4]
            self.snrep=1
        elif space_group_number==79:
            self.name='I4'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[4,2,1]
            self.btype='I'
            self.blt=4
            self.nsg=3
            self.sname=['1','2..','4..']
            self.swycn=[1,1,1]
            self.swyc_mult=[4,2,1]
            self.snrep=3
        elif space_group_number==80:
            self.name='I4_1'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z']
            self.trans=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]]
            self.nwyc=2
            self.wycs=['+x+y+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0]]
            self.wmult=[4,2]
            self.btype='I'
            self.blt=4
            self.nsg=2
            self.sname=['1','2..']
            self.swycn=[1,1]
            self.swyc_mult=[4,2]
            self.snrep=2
        elif space_group_number==81:
            self.name='P-4'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=8
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[4,2,2,2,1,1,1,1]
            self.btype='P'
            self.blt=4
            self.nsg=3
            self.sname=['1','2..','-4..']
            self.swycn=[1,3,4]
            self.swyc_mult=[4,2,1]
            self.snrep=2
        elif space_group_number==82:
            self.name='I-4'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=7
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0.5,0.75],[0,0.5,0.25],[0,0,0.5],[0,0,0]]
            self.wmult=[4,2,2,1,1,1,1]
            self.btype='I'
            self.blt=4
            self.nsg=3
            self.sname=['1','2..','-4..']
            self.swycn=[1,2,4]
            self.swyc_mult=[4,2,1]
            self.snrep=2
        elif space_group_number==83:
            self.name='P4/m'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x-y-z','+x+y-z','+y-x-z','-y+x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=12
            self.wycs=['+x+y+z','+x+y+0','+x+y+0','+0+0+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.5],[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0],
                        [0,0.5,0.5],[0,0.5,0],[0.5,0.5,0.5],[0.5,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,2,2,2,2,1,1,1,1]
            self.btype='P'
            self.blt=4
            self.nsg=6
            self.sname=['1','m..','2..','4..','2/m..','4/m..']
            self.swycn=[1,2,1,2,2,4]
            self.swyc_mult=[8,4,4,2,2,1]
            self.snrep=4
        elif space_group_number==84:
            self.name='P4_2/m'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x-y-z','+x+y-z','+y-x-z','-y+x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5]]
            self.nwyc=11
            self.wycs=['+x+y+z','+x+y+0','+0+0+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0],
                        [0.5,0.5,0.25],[0,0,0.25],[0,0.5,0.5],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=5
            self.sname=['1','m..','2..','-4..','2/m..']
            self.swycn=[1,1,3,2,4]
            self.swyc_mult=[8,4,4,2,2]
            self.snrep=3
        elif space_group_number==85:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4/n_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x-y-z','+x+y-z','+y-x-z','-y+x-z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0]]
                self.nwyc=7
                self.wycs=['+x+y+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.25,0.75,0],[0,0,0.5],[0,0,0],[0.25,0.25,0],[0.25,0.75,0.5],[0.25,0.75,0]]
                self.wmult=[8,4,4,4,2,2,2]
            else:
                self.name='P4/n_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x-y-z','+x+y-z','+y-x-z','-y+x-z']
                self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0]]
                self.nwyc=7
                self.wycs=['+x+y+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0.25,0.25,0.5],[0.25,0.25,0],[0,0.5,0],[0,0,0.5],[0,0,0]]
                self.wmult=[8,4,4,4,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=5
            self.sname=['1','2..','-1','4..','-4..']
            self.swycn=[1,1,2,1,2]
            self.swyc_mult=[8,4,4,2,2]
            self.snrep=4
        elif space_group_number==86:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4_2/n_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x-y-z','+x+y-z','+y-x-z','-y+x-z']
                self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0]]
                self.nwyc=7
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0.5,0],[0.25,0.25,0.75],[0.25,0.25,0.25],[0,0,0.5],[0,0,0]]
                self.wmult=[8,4,4,4,4,2,2]
            else:
                self.name='P_2/n_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x-y-z','+x+y-z','+y-x-z','-y+x-z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5],[0,0,0],[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]]
                self.nwyc=7
                self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.25,0.25,0],[0.75,0.25,0],[0,0,0.5],[0,0,0],[0.25,0.25,0.75],[0.25,0.25,0.25]]
                self.wmult=[8,4,4,4,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=4
            self.sname=['1','2..','-1','-4..']
            self.swycn=[1,2,2,2]
            self.swyc_mult=[8,4,4,2]
            self.snrep=2
        elif space_group_number==87:
            self.name='I4/m'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x-y-z','+x+y-z','+y-x-z','-y+x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=9
            self.wycs=['+x+y+z','+x+y+0','+0+0+z','+0+0+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0],[0.25,0.25,0.25],[0,0,0],[0,0.5,0.25],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,2,2,2,1,1]
            self.btype='I'
            self.blt=4
            self.nsg=8
            self.sname=['1','m..','2..','-1','4..','-4..','2/m..','4/m..']
            self.swycn=[1,1,1,1,1,1,1,2]
            self.swyc_mult=[8,4,4,4,2,2,2,1]
            self.snrep=5
        elif space_group_number==88:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='I4_1/a_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x-y-z','+x+y-z','+y-x-z','-y+x-z']
                self.trans=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75],[0,0.5,0.25],[0.5,0,0.75],[0,0,0],[0.5,0.5,0.5]]
                self.nwyc=6
                self.wycs=['+x+y+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0.25,0.625],[0,0.25,0.125],[0,0,0.5],[0,0,0]]
                self.wmult=[8,4,4,4,2,2]
            else:
                self.name='I4_1/a_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x-y-z','+x+y-z','+y-x-z','-y+x-z']
                self.trans=[[0,0,0],[0.5,0,0.5],[0.75,0.25,0.25],[0.75,0.75,0.75],[0,0,0],[0.5,0,0.5],[0.25,0.75,0.75],[0.25,0.25,0.25]]
                self.nwyc=6
                self.wycs=['+x+y+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.25,0],[0,0,0.5],[0,0,0],[0,0.25,0.625],[0,0.25,0.125]]
                self.wmult=[8,4,4,4,2,2]
            self.btype='I'
            self.blt=4
            self.nsg=4
            self.sname=['1','2..','-1','-4..']
            self.swycn=[1,1,2,2]
            self.swyc_mult=[8,4,4,2]
            self.snrep=2
        elif space_group_number==89:
            self.name='P422'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=16
            self.wycs=['+x+u+z','+x+0+0','+x+0+0','+x+0+0','+x+0+0','+x+x+0','+x+x+0','+0+0+z','+0+0+z','+0+0+z',
                       '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0.5],[0,0.5,0.5],[0,0,0],[0,0,0.5],[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0],
                        [0.5,0,0.5],[0.5,0,0],[0.5,0.5,0.5],[0.5,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,4,4,4,2,2,2,2,1,1,1,1]
            self.btype='P'
            self.blt=3
            self.nsg=7
            self.sname=['1','.2.','..2','2..','4..','222','422']
            self.swycn=[1,4,2,1,2,2,4]
            self.swyc_mult=[8,4,4,4,2,2,1]
            self.snrep=5
        elif space_group_number==90:
            self.name='P42_12'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0]]
            self.nwyc=7
            self.wycs=['+x+y+z','+x+x+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.5],[0,0,0],[0,0,0],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=5
            self.sname=['1','..2','2..','4..','2.2 2']
            self.swycn=[1,2,1,1,2]
            self.swyc_mult=[8,4,4,2,2]
            self.snrep=4
        elif space_group_number==91:
            self.name='P4_122'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0,0.25],[0,0,0.75],[0,0,0],[0,0,0.5],[0,0,0.75],[0,0,0.25]]
            self.nwyc=4
            self.wycs=['+x+y+z','+x+x+0','+0+y+0','+0+y+0']
            self.wtran=[[0,0,0],[0,0,0.375],[0.5,0,0],[0,0,0]]
            self.wmult=[8,4,4,4]
            self.btype='P'
            self.blt=4
            self.nsg=3
            self.sname=['1','..2','.2.']
            self.swycn=[1,1,2]
            self.swyc_mult=[8,4,4]
            self.snrep=3
        elif space_group_number==92:
            self.name='P4_12_12'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0.5],[0.5,0.5,0.25],[0.5,0.5,0.75],[0.5,0.5,0.25],[0.5,0.5,0.75],[0,0,0],[0,0,0.5]]
            self.nwyc=2
            self.wycs=['+x+y+z','+x+x+0']
            self.wtran=[[0,0,0],[0,0,0]]
            self.wmult=[8,4]
            self.btype='P'
            self.blt=4
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,1]
            self.swyc_mult=[8,4]
            self.snrep=2
        elif space_group_number==93:
            self.name='P4_222'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5]]
            self.nwyc=16
            self.wycs=['+x+y+z','+x+x+0','+x+x+0','+x+0+0','+x+0+0','+x+0+0','+x+0+0',
                       '+0+0+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.75],[0,0,0.25],[0,0.5,0],[0,0,0.5],[0,0.5,0.5],[0,0,0],
                        [0,0.5,0],[0.5,0.5,0],[0,0,0],[0.5,0.5,0.25],[0,0,0.25],[0,0.5,0.5],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,4,4,4,4,4,2,2,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=6
            self.sname=['1','..2','.2.','2..','2.2 2','222 .']
            self.swycn=[1,2,4,3,2,4]
            self.swyc_mult=[8,4,4,4,2,2]
            self.snrep=4
        elif space_group_number==94:
            self.name='P4_22_12'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0]]
            self.nwyc=7
            self.wycs=['+x+y+z','+x+x+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0.5],[0,0,0],[0,0.5,0],[0,0,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=4
            self.sname=['1','..2','2..','2.2 2']
            self.swycn=[1,2,2,2]
            self.swyc_mult=[8,4,4,2]
            self.snrep=3
        elif space_group_number==95:
            self.name='P4_322'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0.5],[0,0,0.75],[0,0,0.25],[0,0,0],[0,0,0.5],[0,0,0.25],[0,0,0.75]]
            self.nwyc=4
            self.wycs=['+x+y+z','+x+x+0','+0+y+0','+0+y+0']
            self.wtran=[[0,0,0],[0,0,0.625],[0.5,0,0],[0,0,0]]
            self.wmult=[8,4,4,4]
            self.btype='P'
            self.blt=4
            self.nsg=3
            self.sname=['1','..2','.2.']
            self.swycn=[1,1,2]
            self.swyc_mult=[8,4,4]
            self.snrep=3
        elif space_group_number==96:
            self.name='P4_32_12'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0.5],[0.5,0.5,0.75],[0.5,0.5,0.25],[0.5,0.5,0.75],[0.5,0.5,0.25],[0,0,0],[0,0,0.5]]
            self.nwyc=2
            self.wycs=['+x+y+z','+x+x+0']
            self.wtran=[[0,0,0],[0,0,0]]
            self.wmult=[8,4]
            self.btype='P'
            self.blt=4
            self.nsg=2
            self.sname=['1','..2']
            self.swycn=[1,1]
            self.swyc_mult=[8,4]
            self.snrep=2
        elif space_group_number==97:
            self.name='I422'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=11
            self.wycs=['+x+y+z','+x+x+0','+x+0+0','+x+0+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0.25],[0,0,0.5],[0,0,0],[0,0,0],[0,0.5,0],[0,0,0],[0,0.5,0.25],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,4,2,2,2,1,1]
            self.btype='I'
            self.blt=4
            self.nsg=9
            self.sname=['1','..2','.2.','..2','2..','4..','2.2 2','222 .','422']
            self.swycn=[1,1,2,1,1,1,1,1,2]
            self.swyc_mult=[8,4,4,4,4,2,2,2,1]
            self.snrep=6
        elif space_group_number==98:
            self.name='I4_122'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75],[0.5,0,0.75],[0,0.5,0.25],[0.5,0.5,0.5],[0,0,0]]
            self.nwyc=7
            self.wycs=['+x+y+z','+x+0+0','-x+x+0','+x+x+0','+0+0+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.25,0.125],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2]
            self.btype='I'
            self.blt=4
            self.nsg=5
            self.sname=['1','.2.','..2','2..','2.2 2']
            self.swycn=[1,1,2,1,2]
            self.swyc_mult=[8,4,4,4,2]
            self.snrep=4
        elif space_group_number==99:
            self.name='P4mm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=7
            self.wycs=['+x+y+z','+x+0+z','+x+0+z','+x+x+z','+0+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0,0],[0.5,0,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,4,2,1,1]
            self.btype='P'
            self.blt=4
            self.nsg=5
            self.sname=['1','.m.','..m','2mm','4mm']
            self.swycn=[1,2,1,1,2]
            self.swyc_mult=[8,4,4,2,1]
            self.snrep=5
        elif space_group_number==100:
            self.name='P4bm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=4
            self.wycs=['+x+u+z','+x+x+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0.5,0,0],[0,0,0]]
            self.wmult=[8,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=4
            self.sname=['1','..m','2.m m','4..']
            self.swycn=[1,1,1,1]
            self.swyc_mult=[8,4,2,2]
            self.snrep=4
        elif space_group_number==101:
            self.name="P4_2cm"
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0]]
            self.nwyc=5
            self.wycs=['+x+y+z','+x+x+z','+0+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=4
            self.sname=['1','..m','2..','2.m m']
            self.swycn=[1,1,1,2]
            self.swyc_mult=[8,4,4,2]
            self.snrep=4
        elif space_group_number==102:
            self.name='P4_2nm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0]]
            self.nwyc=4
            self.wycs=['+x+y+z','+x+x+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,2]
            self.btype='P'
            self.blt=4
            self.nsg=4
            self.sname=['1','..m','2..','2.m m']
            self.swycn=[1,1,1,1]
            self.swyc_mult=[8,4,4,2]
            self.snrep=4
        elif space_group_number==103:
            self.name='P4cc'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5]]
            self.nwyc=4
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[8,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=3
            self.sname=['1','2..','4..']
            self.swycn=[1,1,2]
            self.swyc_mult=[8,4,2]
            self.snrep=3
        elif space_group_number==104:
            self.name='P4nc'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[8,4,2]
            self.btype='P'
            self.blt=4
            self.nsg=3
            self.sname=['1','2..','4..']
            self.swycn=[1,1,1]
            self.swyc_mult=[8,4,2]
            self.snrep=3
        elif space_group_number==105:
            self.name='P4_2mc'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5]]
            self.nwyc=6
            self.wycs=['+x+y+z','+x+0+z','+x+0+z','+0+0+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=3
            self.sname=['1','.m.','2mm .']
            self.swycn=[1,2,3]
            self.swyc_mult=[8,4,2]
            self.snrep=3
        elif space_group_number==106:
            self.name='P4_2bc'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0.5],[0.5,0.5,0.5]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[8,4,4]
            self.btype='P'
            self.blt=4
            self.nsg=2
            self.sname=['1','2..']
            self.swycn=[1,2]
            self.swyc_mult=[8,4]
            self.snrep=2
        elif space_group_number==107:
            self.name='I4mm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=5
            self.wycs=['+x+y+z','+x+0+z','+x+x+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,2,1]
            self.btype='I'
            self.blt=4
            self.nsg=5
            self.sname=['1','.m.','..m','2mm','4mm']
            self.swycn=[1,1,1,1,1]
            self.swyc_mult=[8,4,4,2,1]
            self.snrep=5
        elif space_group_number==108:
            self.name='I4cm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5]]
            self.nwyc=4
            self.wycs=['+x+y+z','+x+x+z','+0+0+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0.5,0],[0.5,0,0],[0,0,0]]
            self.wmult=[8,4,2,2]
            self.btype='I'
            self.blt=4
            self.nsg=4
            self.sname=['1','..m','2.m m','4..']
            self.swycn=[1,1,1,1]
            self.swyc_mult=[8,4,2,2]
            self.snrep=4
        elif space_group_number==109:
            self.name='I4_1md'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0.5,0.5,0.5,],[0,0.5,0.25],[0.5,0,0.75],[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75]]
            self.nwyc=3
            self.wycs=['+x+y+z','+0+y+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0]]
            self.wmult=[8,4,2]
            self.btype='I'
            self.blt=4
            self.nsg=3
            self.sname=['1','.m.','2mm .']
            self.swycn=[1,1,1]
            self.swyc_mult=[8,4,2]
            self.snrep=3
        elif space_group_number==110:
            self.name='I4_1cd'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75],[0,0,0.5],[0.5,0.5,0],[0,0.5,0.75],[0.5,0,0.25]]
            self.nwyc=2
            self.wycs=['+x+y+z','+0+0+z']
            self.wtran=[[0,0,0],[0,0,0]]
            self.wmult=[8,4]
            self.btype='I'
            self.blt=4
            self.nsg=2
            self.sname=['1','2..']
            self.swycn=[1,1]
            self.swyc_mult=[8,4]
            self.snrep=2
        elif space_group_number==111:
            self.name='P-42m'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','-x+y-z','+x-y-z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=15
            self.wycs=['+x+y+z','+x+x+z','+0+0+z','+x+0+0','+x+0+0','+x+0+0','+x+0+0','+0+0+z','+0+0+z',
                       '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0],[0,0.5,0],[0,0,0.5],[0,0.5,0.5],[0,0,0],[0.5,0.5,0],[0,0,0],
                        [0.5,0,0.5],[0.5,0,0],[0.5,0.5,0],[0,0,0.5],[0.5,0.5,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,4,4,2,2,2,2,1,1,1,1]
            self.btype='P'
            self.blt=4
            self.nsg=7
            self.sname=['1','..m','2..','.2.','2.m m','222 .','-42m']
            self.swycn=[1,1,1,4,2,2,4]
            self.swyc_mult=[8,4,4,4,2,2,1]
            self.snrep=5
        elif space_group_number==112:
            self.name='P-42c'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','-x+y-z','+x-y-z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5]]
            self.nwyc=14
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z','+0+y+0','+x+0+0','+0+y+0','+x+0+0',
                       '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0.25],[0,0.5,0.25],[0.5,0,0.25],[0,0,0.25],
                        [0.5,0.5,0],[0,0,0],[0,0.5,0.25],[0.5,0.5,0.25],[0.5,0,0.25],[0,0,0.25]]
            self.wmult=[8,4,4,4,4,4,4,4,2,2,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=5
            self.sname=['1','2..','.2.','-4..','222 .']
            self.swycn=[1,3,4,2,4]
            self.swyc_mult=[8,4,4,2,2]
            self.snrep=3
        elif space_group_number==113:
            self.name='P-42_1m'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','-x+y-z','+x-y-z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=6
            self.wycs=['+x+y+z','+x+x+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=5
            self.sname=['1','..m','2..','2.m m','-4..']
            self.swycn=[1,1,1,1,2]
            self.swyc_mult=[8,4,4,2,2]
            self.snrep=4
        elif space_group_number==114:
            self.name='P-42_1c'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','-x+y-z','+x-y-z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]
            self.nwyc=5
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=3
            self.sname=['1','2..','-4..']
            self.swycn=[1,2,2]
            self.swyc_mult=[8,4,2]
            self.snrep=2
        elif space_group_number==115:
            self.name='P-4m2'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=12
            self.wycs=['+x+y+z','+x+0+z','+x+0+z','+x+x+0','+x+x+0','+0+0+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0,0.5],[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0.5],[0.5,0.5,0.5],[0.5,0.5,0],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2,1,1,1,1]
            self.btype='P'
            self.blt=4
            self.nsg=5
            self.sname=['1','.m.','..2','2mm .','-4m2']
            self.swycn=[1,2,2,3,4]
            self.swyc_mult=[8,4,4,2,1]
            self.snrep=4
        elif space_group_number==116:
            self.name='P-4c2'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5]]
            self.nwyc=10
            self.wycs=['+x+y+z','+0+0+z','+0+0+z','+0+0+z','+x+x+0','+x+x+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0.75],[0,0,0.25],[0.5,0.5,0],[0,0,0],[0.5,0.5,0.25],[0,0,0.25]]
            self.wmult=[8,4,4,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=5
            self.sname=['1','2..','..2','-4..','2.2 2']
            self.swycn=[1,3,2,2,2]
            self.swyc_mult=[8,4,4,2,2]
            self.snrep=3
        elif space_group_number==117:
            self.name='P-4b2'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=9
            self.wycs=['+x+y+z','+x+x+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0.5],[0,0.5,0],[0,0.5,0],[0,0,0],[0,0.5,0.5],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=5
            self.sname=['1','..2','2..','2.2 2','-4..']
            self.swycn=[1,2,2,2,2]
            self.swyc_mult=[8,4,4,2,2]
            self.snrep=3
        elif space_group_number==118:
            self.name='P-4n2'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]
            self.nwyc=9
            self.wycs=['+x+y+z','+0+0+z','+x+x+0','+x+x+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0.5,0.25],[0,0.5,0.25],[0,0,0],[0,0.5,0.75],[0,0.5,0.25],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=6
            self.sname=['1','2..','..2','2..','2.2 2','-4..']
            self.swycn=[1,1,2,1,2,2]
            self.swyc_mult=[8,4,4,4,2,2]
            self.snrep=4
        elif space_group_number==119:
            self.name='I-4m2'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=10
            self.wycs=['+x+y+z','+x+0+z','+x+x+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0.25],[0,0,0],[0,0.5,0],[0,0,0],[0,0.5,0.75],[0,0.5,0.25],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,2,2,1,1,1,1]
            self.btype='I'
            self.blt=4
            self.nsg=5
            self.sname=['1','.m.','..2','2mm .','-4m2']
            self.swycn=[1,1,2,2,4]
            self.swyc_mult=[8,4,4,2,1]
            self.snrep=4
        elif space_group_number==120:
            self.name='I-4c2'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','+y+x-z','-y-x-z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5]]
            self.nwyc=9
            self.wycs=['+x+y+z','+x+x+0','+0+0+z','+0+0+z','+x+x+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0.5,0],[0,0,0],[0,0,0.25],[0,0.5,0],[0,0.5,0.25],[0,0,0],[0,0,0.25]]
            self.wmult=[8,4,4,4,4,2,2,2,2]
            self.btype='I'
            self.blt=4
            self.nsg=7
            self.sname=['1','..2','2..','..2','2.2 2','-4..' '2.2 2']
            self.swycn=[1,1,2,1,1,2,1]
            self.swyc_mult=[8,4,4,4,2,2,2]
            self.snrep=4
        elif space_group_number==121:
            self.name='I-42m'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','-x+y-z','+x-y-z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=10
            self.wycs=['+x+y+z','+x+x+z','+0+0+z','+x+0+0','+x+0+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0],[0,0,0.5],[0,0,0],[0,0,0],[0,0.5,0.25],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,4,4,2,2,2,1,1]
            self.btype='I'
            self.blt=4
            self.nsg=8
            self.sname=['1','..m','2..','.2.','2.m m','-4..','222 .','-42m']
            self.swycn=[1,1,1,2,1,1,1,2]
            self.swyc_mult=[8,4,4,4,2,2,2,1]
            self.snrep=5
        elif space_group_number==122:
            self.name='I-42d'
            self.pgsym=['+x+y+z','-x-y+z','+y-x-z','-y+x-z','-x+y-z','+x-y-z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0,0.75],[0.5,0,0.75],[0.5,0,0.75],[0.5,0,0.75]]
            self.nwyc=5
            self.wycs=['+x+y+z','+x+0+0','+0+0+z','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.25,0.125],[0,0,0],[0,0,0.5],[0,0,0]]
            self.wmult=[8,4,4,2,2]
            self.btype='I'
            self.blt=4
            self.nsg=4
            self.sname=['1','.2.','2..','-4..']
            self.swycn=[1,1,1,2]
            self.swyc_mult=[8,4,4,2]
            self.snrep=3
        elif space_group_number==123:
            self.name='P4/mmm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],
                        [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=21
            self.wycs=['+x+y+z','+x+0+z','+x+0+z','+x+x+z','+x+y+0','+x+y+0',
                       '+x+0+0','+x+0+0','+x+0+0','+x+0+0','+x+x+0','+x+x+0',
                       '+0+0+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0],
                        [0,0.5,0.5],[0,0.5,0],[0,0,0.5],[0,0,0],[0,0,0.5],[0,0,0],
                        [0,0.5,0],[0.5,0.5,0],[0,0,0],[0,0.5,0],[0,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[16,8,8,8,8,8,4,4,4,4,4,4,4,2,2,2,2,1,1,1,1]
            self.btype='P'
            self.blt=4
            self.nsg=10
            self.sname=['1','.m.','..m','m..','m2m .','m.2 m','2mm .','4mm','mmm .','4/mmm']
            self.swycn=[1,2,1,2,4,2,1,2,2,4]
            self.swyc_mult=[16,8,8,8,4,4,4,2,2,1]
            self.snrep=8
        elif space_group_number==124:
            self.name='P4/mcc'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],
                        [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5]]
            self.nwyc=14
            self.wycs=['+x+y+z','+x+y+0','+x+0+0','+x+0+0','+x+x+0',
                       '+0+0+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0.25],[0,0,0.25],[0,0,0.25],
                        [0,0.5,0],[0.5,0.5,0],[0,0,0],[0,0.5,0.25],[0,0.5,0],[0.5,0.5,0],[0.5,0.5,0.25],[0,0,0],[0,0,0.25]]
            self.wmult=[16,8,8,8,8,8,4,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=12
            self.sname=['1','m..','.2.','..2','2..','4..','222 .','2/m..','4/m..','422','4/m..','422']
            self.swycn=[1,1,2,1,1,2,1,1,1,1,1,1]
            self.swyc_mult=[16,8,8,8,8,4,4,4,2,2,2,2]
            self.snrep=6
        elif space_group_number==125:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4/nbm_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],
                            [0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0]]
                self.nwyc=14
                self.wycs=['+x+y+z','+x+x+z','+x+0+0','+x+0+0','+x+x+0','+x+x+0','+0+0+z','+0+0+z',
                           '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.5,0],[0,0,0.5],[0,0,0],[0,0,0.5],[0,0,0],[0,0.5,0],[0,0,0],
                            [0.25,0.25,0.5],[0.25,0.25,0],[0,0.5,0.5],[0,0.5,0],[0,0,0.5],[0,0,0]]
                self.wmult=[16,8,8,8,8,8,4,4,4,4,2,2,2,2]
            else:
                self.name='P4/nbm_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0],[0.5,0.5,0],
                            [0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0],[0.5,0.5,0]]
                self.nwyc=14
                self.wycs=['+x+y+z','+x-x+z','+x+0+0','+x+0+0','+x+x+0','+x+x+0','+0+0+z','+0+0+z',
                           '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0.25,0.5],[0,0.25,0],[0,0,0.5],[0,0,0],[0.75,0.25,0],[0.25,0.25,0],
                            [0,0,0.5],[0,0,0],[0.75,0.25,0.5],[0.75,0.25,0],[0.25,0.25,0.5],[0.25,0.25,0]]
                self.wmult=[16,8,8,8,8,8,4,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=9
            self.sname=['1','..m','.2.','..2','2.m m','4..','..2/m','-42m','422']
            self.swycn=[1,1,2,2,1,1,2,2,2]
            self.swyc_mult=[16,8,8,8,4,4,4,2,2]
            self.snrep=6
        elif space_group_number==126:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4/nnc_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],
                            [0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]
                self.nwyc=11
                self.wycs=['+x+y+z','+x+0+0','+x+0+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0.5],[0,0,0],[0,0,0],[0.5,0,0],[0,0,0],[0.25,0.25,0.25],[0.5,0,0.25],[0.5,0,0],[0,0,0.5],[0,0,0]]
                self.wmult=[16,8,8,8,8,4,8,4,4,2,2]
            else:
                self.name='P4/nnc_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0,0,0.5],[0.5,0.5,0.5],
                            [0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0,0,0.5],[0.5,0.5,0.5]]
                self.nwyc=11
                self.wycs=['+x+y+z','+x+0+0','+x+0+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.75,0.25],[0,0.25,0.25],[0,0,0.25],[0.25,0.75,0],[0.25,0.25,0],[0,0,0],[0.25,0.75,0],[0.25,0.75,0.75],[0.25,0.25,0.75],[0.25,0.25,0.25]]
                self.wmult=[16,8,8,8,8,4,8,4,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=9
            self.sname=['1','.2.','..2','2..','4..','=1','-4..','222 .','422']
            self.swycn=[1,2,1,1,1,1,1,1,2]
            self.swyc_mult=[16,8,8,8,4,8,4,4,2]
            self.snrep=5
        elif space_group_number==127:
            self.name='P4/mbm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],
                        [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0]]
            self.nwyc=12
            self.wycs=['+x+y+z','+x+x+z','+x+y+0','+x+y+0','+x+x+0','+x+x+0',
                       '+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0.5],[0,0,0],[0,0.5,0.5],[0,0.5,0],
                        [0,0.5,0],[0,0,0],[0,0.5,0],[0,0.5,0.5],[0,0,0.5],[0,0,0]]
            self.wmult=[16,8,8,8,4,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=8
            self.sname=['1','..m','m..','m.2 m','2.m m','4..','m.m m','4/m..']
            self.swycn=[1,1,2,2,1,1,2,2]
            self.swyc_mult=[16,8,8,4,4,4,2,2]
            self.snrep=6
        elif space_group_number==128:
            self.name='P4/mnc'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],
                        [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]
            self.nwyc=9
            self.wycs=['+x+y+z','+x+y+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0.25],[0,0.5,0],[0,0,0],[0,0.5,0.25],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[16,8,8,8,4,4,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=8
            self.sname=['1','m..','..2','2..','4..','2.2 2','2/m..','4/m..']
            self.swycn=[1,1,1,1,1,1,1,2]
            self.swyc_mult=[16,8,8,8,4,4,4,2]
            self.snrep=5
        elif space_group_number==129:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4/nmm_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0],
                            [0.5,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0]]
                self.nwyc=11
                self.wycs=['+x+y+z','+x+x+z','+0+y+z','+x+x+0','+x+x+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0,0.5],[0,0,0],[0,0,0],[0.25,0.25,0.5],[0.25,0.25,0],[0,0.5,0],[0,0,0.5],[0,0,0]]
                self.wmult=[16,8,8,8,8,4,4,4,2,2,2]
            else:
                self.name='P4/nmm_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0.5,0],[0.5,0,0],[0.5,0.5,0],[0,0,0],
                            [0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0.5,0],[0.5,0,0],[0.5,0.5,0],[0,0,0]]
                self.nwyc=11
                self.wycs=['+x+y+z','+x+x+z','+0+y+z','+x+x+0','+x+x+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0.25,0,0],[0,0,0.5],[0,0,0],[0.75,0.25,0],[0,0,0.5],[0,0,0],[0.25,0.25,0],[0.75,0.25,0.5],[0.75,0.25,0]]
                self.wmult=[16,8,8,8,8,4,4,4,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=8
            self.sname=['1','..m','.m.','..2','2mm .','..2/m','4mm','-4m2']
            self.swycn=[1,1,1,2,1,2,1,2]
            self.swyc_mult=[16,8,8,8,4,4,2,2]
            self.snrep=5
        elif space_group_number==130:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4/ncc_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0.5],[0,0,0.5],
                            [0.5,0.5,0],[0.5,0.5,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]
                self.nwyc=7
                self.wycs=['+x+y+z','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0.25],[0,0,0],[0,0.5,0],[0.25,0.25,0],[0,0,0],[0,0,0.25]]
                self.wmult=[16,8,8,4,8,4,4]
            else:
                self.name='P4/ncc_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0.5],[0,0,0.5],
                            [0,0,0],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0.5],[0,0,0.5]]
                self.nwyc=7
                self.wycs=['+x+y+z','+x-x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0.25],[0.75,0.25,0],[0.25,0.25,0],[0,0,0],[0.75,0.25,0],[0.75,0.25,0.25]]
                self.wmult=[16,8,8,4,8,4,4]
            self.btype='P'
            self.blt=4
            self.nsg=7
            self.sname=['1','..2','2..','4..','=1','-4..','2.2 2']
            self.swycn=[1,1,1,1,1,1,1]
            self.swyc_mult=[16,8,8,4,8,4,4]
            self.snrep=4
        elif space_group_number==131:
            self.name='P4_2/mmc'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5]]
            self.nwyc=18
            self.wycs=['+x+y+z','+x+y+0','+0+y+z','+0+y+z','+x+x+0','+x+0+0','+x+0+0','+x+0+0','+x+0+0','+0+0+z','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0.5,0,0],[0,0,0],[0,0,0.25],[0,0.5,0],[0,0,0.5],[0,0.5,0.5],[0,0,0],[0,0.5,0],[0.5,0.5,0],[0,0,0],
                        [0.5,0.5,0.25],[0,0,0.25],[0,0.5,0.5],[0,0.5,0],[0.5,0.5,0],[0,0,0]]
            self.wmult=[16,8,8,8,8,4,4,4,4,4,4,4,2,2,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=8
            self.sname=['1','m..','.m.','..2','m2m .','2mm .','-4m2','mmm .']
            self.swycn=[1,1,2,1,4,3,2,4]
            self.swyc_mult=[16,8,8,8,4,4,2,2]
            self.snrep=6
        elif space_group_number==132:
            self.name='P4_2/mcm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0],
                        [0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0],[0,0,0]]
            self.nwyc=16
            self.wycs=['+x+y+z','+x+x+z','+x+y+0','+x+0+0','+x+0+0','+0+0+z','+x+x+0','+x+x+0',
                       '+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0.5,0.25],[0,0,0.25],[0,0.5,0],[0,0,0.5],[0,0,0],
                        [0.5,0.5,0],[0,0,0],[0,0.5,0],[0,0.5,0.25],[0.5,0.5,0.25],[0.5,0.5,0],[0,0,0.25],[0,0,0]]
            self.wmult=[16,8,8,8,8,8,4,4,4,4,4,4,2,2,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=13
            self.sname=['1','..m','m..','.2.','2..','m.2 m','2.m m','2/m..','222 .','-42m','m.m m','-4m','m.m m']
            self.swycn=[1,1,1,2,1,2,2,1,1,1,1,1,1]
            self.swyc_mult=[16,8,8,8,8,4,4,4,4,2,2,2,2]
            self.snrep=7
        elif space_group_number==133:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4_2/nbc_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0.5],[0,0,0.5],[0.5,0.5,0],[0.5,0.5,0],
                            [0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0],[0.5,0.5,0],[0.5,0.5,0],[0,0,0.5],[0,0,0.5]]
                self.nwyc=11
                self.wycs=['+x+y+z','+x+x+0','+x+0+0','+x+0+0','+0+0+z','+0+0+z',
                           '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.5,0],[0,0,0.75],[0,0,0.25],[0,0,0],[0,0.5,0],
                            [0.25,0.25,0.25],[0,0,0],[0,0.5,0],[0,0,0.25],[0,0.5,0.25]]
                self.wmult=[16,8,8,8,8,8,8,4,4,4,4]
            else:
                self.name='P4_2/nbc_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0,0],[0,0.5,0],[0,0,0.5],[0.5,0.5,0.5],
                            [0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0,0],[0,0.5,0],[0,0,0.5],[0.5,0.5,0.5]]
                self.nwyc=11
                self.wycs=['+x+y+z','+x+x+0','+x+0+0','+x+0+0','+0+0+z','+0+0+z',
                           '+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0.25],[0,0.25,0.5],[0,0.25,0],[0.75,0.25,0],[0.25,0.25,0],
                            [0,0,0],[0.75,0.25,0.75],[0.25,0.25,0.25],[0.75,0.25,0],[0.25,0.25,0]]
                self.wmult=[16,8,8,8,8,8,8,4,4,4,4]
            self.btype='P'
            self.blt=4
            self.nsg=8
            self.sname=['1','..2','.2.','2..','-1','-4..','2.2 2','222 .']
            self.swycn=[1,1,2,2,1,1,1,2]
            self.swyc_mult=[16,8,8,8,8,4,4,4]
            self.snrep=4
        elif space_group_number==134:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4_2/nnm_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],
                            [0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0]]
                self.nwyc=14
                self.wycs=['+x+y+z','+x+x+z','+x+x+0','+x+x+0','+x+0+0','+x+0+0',
                           '+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0.5,0.75],[0,0.5,0.25],[0,0,0.5],[0,0,0],
                            [0,0.5,0],[0,0,0],[0.75,0.75,0.75],[0.25,0.25,0.25],[0,0.5,0.25],[0,0.5,0],[0,0,0.5],[0,0,0]]
                self.wmult=[16,8,8,8,8,8,8,4,4,4,4,4,2,2]
            else:
                self.name='P4_2/nnm_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0,0.5],[0,0.5,0.5],[0,0,0],[0.5,0.5,0],
                            [0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0,0.5],[0,0.5,0.5],[0,0,0],[0.5,0.5,0]]
                self.nwyc=14
                self.wycs=['+x+y+z','+x-x+z','+x+x+0','+x+x+0','+x+0+0','+x+0+0',
                           '+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0],[0,0.25,0.25],[0,0.25,0.75],
                            [0.25,0.25,0],[0.75,0.25,0],[0,0,0],[0,0,0.5],[0.25,0.25,0],[0.25,0.25,0.25],[0.75,0.25,0.25],[0.25,0.75,0.25]]
                self.wmult=[16,8,8,8,8,8,8,4,4,4,4,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=10
            self.sname=['1','..m','..2','.2.','2..','2.m m','..2/m','2.2 2','222 .','-42m']
            self.swycn=[1,1,2,2,1,1,2,1,1,2]
            self.swyc_mult=[16,8,8,8,8,4,4,4,4,2]
            self.snrep=6
        elif space_group_number==135:
            self.name='P4_2/mbc'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0.5],[0.5,0.5,0.5],
                        [0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0.5,0.5,0],[0.5,0.5,0],[0.5,0.5,0.5],[0.5,0.5,0.5]]
            self.nwyc=9
            self.wycs=['+x+y+z','+x+y+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0.5,0.25],[0,0.5,0],[0,0,0],[0,0.5,0.25],[0,0.5,0],[0,0,0.25],[0,0,0]]
            self.wmult=[16,8,8,8,8,4,4,4,4]
            self.btype='P'
            self.blt=4
            self.nsg=8
            self.sname=['1','m..','..2','2..','2.2 2','2/m..','-4..','2/m..']
            self.swycn=[1,1,1,2,1,1,1,1]
            self.swyc_mult=[16,8,8,8,4,4,4,4]
            self.snrep=4
        elif space_group_number==136:
            self.name='P4_2/mnm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0],
                        [0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0]]
            self.nwyc=11
            self.wycs=['+x+y+z','+x+x+z','+x+y+0','+0+0+z','+x-x+0','+x+x+0',
                       '+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0.5,0],[0,0,0],[0,0,0],
                        [0,0,0],[0,0.5,0.25],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[16,8,8,8,4,4,4,4,4,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=9
            self.sname=['1','..m','m..','2..','m.2 m','2.m m','-4..','2/m..','m.m m']
            self.swycn=[1,1,1,1,2,1,1,1,2]
            self.swyc_mult=[16,8,8,8,4,4,4,4,2]
            self.snrep=6
        elif space_group_number==137:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4_2/nmc_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0],
                            [0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5]]
                self.nwyc=8
                self.wycs=['+x+y+z','+0+y+z','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0.5,0],[0,0,0],[0.25,0.25,0.25],[0,0,0.5],[0,0,0]]
                self.wmult=[16,8,8,4,4,8,2,2]
            else:
                self.name='P4_2/nmc_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0,0.5,0],[0.5,0,0],[0.5,0.5,0.5],[0,0,0.5],
                            [0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0,0.5,0],[0.5,0,0],[0.5,0.5,0.5],[0,0,0.5]]
                self.nwyc=8
                self.wycs=['+x+y+z','+0+y+z','+x-x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0.25,0,0],[0,0,0.25],[0.25,0.25,0],[0.75,0.25,0],[0,0,0],[0.75,0.25,0.25],[0.75,0.25,0.75]]
                self.wmult=[16,8,8,4,4,8,2,2]
            self.btype='P'
            self.blt=4
            self.nsg=6
            self.sname=['1','.m.','..2','2mm .','-1','-4m2']
            self.swycn=[1,1,1,2,1,2]
            self.swyc_mult=[16,8,8,4,8,2]
            self.snrep=4
        elif space_group_number==138:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='P4_2/ncm_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0],[0.5,0.5,0],[0,0,0.5],[0,0,0.5],
                            [0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0.5,0.5,0],[0.5,0.5,0]]
                self.nwyc=10
                self.wycs=['+x+y+z','+x+x+z','+x+x+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.5,0],[0,0,0.75],[0,0,0.25],[0,0,0],[0,0.5,0],[0.25,0.25,0.75],[0.25,0.25,0.25],[0,0,0],[0,0,0.25]]
                self.wmult=[16,8,8,8,8,4,4,4,4,4]
            else:
                self.name='P4_2/ncm_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0,0,0],
                            [0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0,0,0]]
                self.nwyc=10
                self.wycs=['+x+y+z','+x+x+z','+x+x+0','+x+x+0','+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0.75,0.25,0],[0.25,0.25,0],[0,0,0],[0,0,0.5],[0.75,0.25,0.75],[0.75,0.25,0]]
                self.wmult=[16,8,8,8,8,4,4,4,4,4]
            self.btype='P'
            self.blt=4
            self.nsg=8
            self.sname=['1','..m','..2','2..','2.m m','..2/m','-4..','2.2 2']
            self.swycn=[1,1,2,1,1,2,1,1]
            self.swyc_mult=[16,8,8,8,4,4,4,4]
            self.snrep=5
        elif space_group_number==139:
            self.name='I4/mmm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],
                        [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
            self.nwyc=15
            self.wycs=['+x+y+z','+0+y=z','+x+x+z','+x+y+0','+x+x+0','+x+0+0','+x+0+0','+x+x+0',
                       '+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0.5,0.25],[0,0.5,0],[0,0,0],[0,0,0],
                        [0,0.5,0],[0,0,0],[0.25,0.25,0.25],[0,0.5,0.25],[0,0.5,0],[0,0,0.5],[0,0,0]]
            self.wmult=[16,8,8,8,8,4,4,4,4,2,4,2,2,1,1]
            self.btype='I'
            self.blt=4
            self.nsg=13
            self.sname=['1','.m.','..m','m..','..2','m2m .','m.2 m','2mm .','4mm','..2/m','-4m2','mmm .','4/mmm']
            self.swycn=[1,1,1,1,1,2,1,1,1,1,1,1,2]
            self.swyc_mult=[16,8,8,8,8,4,4,4,2,4,2,2,1]
            self.snrep=9
        elif space_group_number==140:
            self.name='I4/mcm'
            self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                        '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
            self.trans=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5],
                        [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0.5],[0,0,0.5],[0,0,0.5],[0,0,0.5]]
            self.nwyc=13
            self.wycs=['+x+y+z','+x+x+z','+x+y+0','+x+0+0','+x+x+0','+x+x+0',
                       '+0+0+z','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
            self.wtran=[[0,0,0],[0,0.5,0],[0,0,0],[0,0,0.25],[0,0,0.25],[0,0.5,0],
                        [0,0.5,0],[0,0,0],[0.25,0.25,0.25],[0,0.5,0],[0,0,0],[0,0.5,0.25],[0,0,0.25]]
            self.wmult=[16,8,8,8,8,4,4,4,4,2,2,2,2]
            self.btype='I'
            self.blt=4
            self.nsg=13
            self.sname=['1','..m','m..','.2.','..2','m.2 m','2.m m','4..','..2/m','m.m m','4/m..','-42m','422']
            self.swycn=[1,1,1,1,1,1,1,1,1,1,1,1,1]
            self.swyc_mult=[16,8,8,8,8,4,4,4,4,2,2,2,2]
            self.snrep=8
        elif space_group_number==141:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='I4_1/amd_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75],[0.5,0,0.75],[0,0.5,0.25],[0.5,0.5,0.5],[0,0,0],
                            [0,0.5,0.25],[0.5,0,0.75],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0.5],[0,0,0],[0.5,0,0.75],[0,0.5,0.25]]
                self.nwyc=9
                self.wycs=['+x+y+z','+0+y+z','+x+x+0','+x+0+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0,0],[0,0.25,0.125],[0,0,0],[0,0.25,0.625],[0,0.25,0.125],[0,0,0.5],[0,0,0]]
                self.wmult=[16,8,8,8,4,4,4,2,2]
            else:
                self.name='I4_1/amd_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0,0.5],[0.25,0.75,0.25],[0.25,0.25,0.75],[0.5,0,0.5],[0,0,0],[0.25,0.75,0.25],[0.25,0.25,0.75],
                            [0,0,0],[0.5,0,0.5],[0.75,0.25,0.75],[0.75,0.75,0.25],[0.5,0,0.5],[0,0,0],[0.75,0.25,0.75],[0.75,0.75,0.25]]
                self.nwyc=9
                self.wycs=['+x+y+z','+0+y+z','+x+x+0','+x+0+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0],[0,0.25,0.875],[0,0,0],[0,0.25,0],[0,0,0.5],[0,0,0],[0,0.25,0.375],[0,0.75,0.125]]
                self.wmult=[16,8,8,8,4,4,4,2,2]
            self.btype='I'
            self.blt=4
            self.nsg=7
            self.sname=['1','.m.','..2','.2.','2mm .','.2/m.','-4m2']
            self.swycn=[1,1,1,1,1,2,2]
            self.swyc_mult=[16,8,8,8,4,4,2]
            self.snrep=5
        elif space_group_number==142:
            l=setting[alternative_setting]
            if l<0.5:         
                self.name='I4_1/acd_oc_1'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0.5,0.5],[0,0.5,0.25],[0.5,0,0.75],[0.5,0,0.25],[0,0.5,0.75],[0.5,0.5,0],[0,0,0.5],
                            [0,0.5,0.25],[0.5,0,0.75],[0,0,0],[0.5,0.5,0.5],[0.5,0.5,0],[0,0,0.5],[0.5,0,0.25],[0,0.5,0.75]]
                self.nwyc=7
                self.wycs=['+x+y+z','+x+x+0','+0+y+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0,0.25],[0.25,0,0.125],[0,0,0],[0,0.25,0.125],[0,0,0.25],[0,0,0]]
                self.wmult=[16,8,8,8,8,4,4]
            else:
                self.name='I4_1/acd_oc_2'
                self.pgsym=['+x+y+z','-x-y+z','-y+x+z','+y-x+z','-x+y-z','+x-y-z','+y+x-z','-y-x-z',
                            '-x-y-z','+x+y-z','+y-x-z','-y+x-z','+x-y+z','-x+y+z','-y-x+z','+y+x+z']
                self.trans=[[0,0,0],[0.5,0,0.5],[0.25,0.75,0.25],[0.25,0.25,0.75],[0.5,0,0],[0,0,0.5],[0.25,0.75,0.75],[0.25,0.25,0.25],
                            [0,0,0],[0.5,0,0.5],[0.75,0.25,0.75],[0.75,0.75,0.25],[0.5,0,0],[0,0,0.5],[0.75,0.25,0.25],[0.75,0.75,0.75]]
                self.nwyc=7
                self.wycs=['+x+y+z','+x+x+0','+x+0+0','+0+0+z','+0+0+0','+0+0+0','+0+0+0']
                self.wtran=[[0,0,0],[0,0.25,0.125],[0,0,0.25],[0,0.25,0],[0,0,0],[0,0.25,0.125],[0,0.25,0.375]]
                self.wmult=[16,8,8,8,8,4,4]
            self.btype='I'
            self.blt=4
            self.nsg=7
            self.sname=['1','..2','.2.','2..','-1','2.2 2','-4..']
            self.swycn=[1,1,1,1,1,1,1]
            self.swyc_mult=[16,8,8,8,8,4,4]
            self.snrep=4
        else:
            self.type=-1
        
        self.op=[] #Storing of the rotation of the molecule in matrix format
        for i in range (0,len(self.pgsym)):
            ns=symop(self.pgsym[i])
            self.op.append(ns.mat)    
        self.wop=[]
        for i in range (0,len(self.wycs)):
            ns=symop(self.wycs[i])
            self.wop.append(ns.mat)
        if self.btype=='P':
            self.btran=[[0,0,0]]
            self.bmult=1
        if self.btype=='I':
            self.btran=[[0,0,0],[0.5,0.5,0.5]]
            self.bmult=2
        if self.btype=='C':
            self.btran=[[0,0,0],[0.5,0.5,0]]
            self.bmult=2
        if self.btype=='F':
            self.btran=[[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]]
            self.bmult=4
        if self.btype=='A':
            self.btran=[[0,0,0],[0,0.5,0.5]]
            self.bmult=2
        self.mlist=[]
    def wycgen(self,wn):
        '''
        Produces a position of wyckoff position no. wn
        '''
        v=[random.random(),random.random(),random.random()]
        nxyz=[0,0,0]
        for i in range (0,3):
            for j in range (0,3):
                nxyz[i]+=self.wop[wn][i][j]*v[j]
        for i in range (0,3):
            nxyz[i]+=self.wtran[wn][i]
        return nxyz
    def wycarr(self,nmpc):
        '''
        Produces a list of Wyckoff position arrangements that add up to nmpc of molecules per cell
        Need to update wycarr!!!
        '''
        mlist=[]
        for i in range (0,nmpc+1):
            mlist.append([]) #Initializing mlist
        mlist[0]=[[]] #Provides the nil solution
        for i in range (0,self.nwyc):
            for j in range (0,nmpc+1-self.wmult[i]*self.bmult):
                for k in range (0,len(mlist[j])):
                    mlist[j+self.wmult[i]*self.bmult].append(mlist[j][k]+[i])
        return mlist[nmpc]

    def wyckoff_preparation (self,nmpc):
        '''
        Prepare a list of Wyckoff position for placing the molecule
        '''
        self.mlist=[]
        for i in range (0,nmpc+1):
            self.mlist.append([])
        self.mlist[0]=[[]]
        for i in range (0,self.nsg):
            if self.snrep>i:
            #A full knapsack for the repeating Wyckoff Positions
                for j in range (0,nmpc+1-self.swyc_mult[i]*self.bmult):
                    for k in range (0,len(self.mlist[j])):
                        self.mlist[j+self.swyc_mult[i]*self.bmult].append(self.mlist[j][k]+[i])
            else:
            #A limited knapsack for the non-repeating Wyckoff positions
                for j in range (nmpc-self.swyc_mult[i]*self.bmult,-1,-1):
                        #                    print "hello j=",j, mlist[j]
                    for l in range (1,1+min(self.swycn[i],(nmpc-j)/(self.swyc_mult[i]*self.bmult))):
                        for k in range (0,len(self.mlist[j])):
                            self.mlist[j+l*self.swyc_mult[i]*self.bmult].append(self.mlist[j][k]+[i for ll in range (0,l)])
#        print self.mlist[nmpc]
        return not (len(self.mlist[nmpc])==0)

    def wyckoff_selection (self,nmpc):
        ll=int(random.uniform(0,len(self.mlist[nmpc])))
        slist=self.mlist[nmpc][ll]
        wlist=[]
        for i in range (0,len(slist)):
            if slist[i]<self.snrep:
                wn=int(random.uniform(0,self.swycn[slist[i]]))
                for j in range (0,slist[i]):
                    wn+=self.swycn[j]
                wlist.append(wn)
            else:
                while True:
                    wn=int(random.uniform(0,self.swycn[slist[i]]))
                    for j in range (0,slist[i]):
                        wn+=self.swycn[j]
                    if not (wn in wlist):
                        break
                wlist.append(wn)
        return wlist

    def wyckoff_counter (self,wyckoff_list):
        '''
        Count the sum of molecules required by the wyckoff list specified
        '''
        counter = 0
        for wyc in wyckoff_list:
            if wyc >= len(self.wmult): #Unsuitable list
                return -1
            counter += self.wmult[wyc]*self.bmult
        return counter

    def get_bravais_system_type(self):
        return self.blt

def allowed_sg_wyckoff_list(nmpc,wyckoff_list):
    '''
    Given a list of Wyckoff position number
    Returns a list of space group that gives the desirable nmpc
    '''
    result = []
    for i in range (0,MAX_AVAILABLE_SPACE_GROUP):
        sg = Sgroup(i)
        if nmpc==sg.wyckoff_counter(wyckoff_list):
            result.append(i)
    return result

def allowed_sg_nmpc(nmpc):
    '''
    Given nmpc, returns the list of space group that can yield the number
    '''
    result = []
    for i in range (0,MAX_AVAILABLE_SPACE_GROUP):
        sg = Sgroup(i)
        if sg.wyckoff_preparation(nmpc):
            result.append(i)
    return result

def select_chiral_sg (sg_list):
    '''
    Select and return the chiral space groups in the given sg_list
    sg_list should be a list of integers
    list of chiral space groups: 
    http://www-chimie.u-strasbg.fr/csd.doc/ConQuest/PortableHTML/conquest_portable-3-325.html
    '''
    result = []
    for sg in sg_list:
        if sg in [1]+range(3,6)+range(16,25)+range(75,81)+range(89,99)\
            +range(143,147)+range(149,156)+range(168,174)+range(177,183)\
            +range(195,200)+range(207,215):
            result.append(sg)
    return result

def select_racemic_sg (sg_list):

    '''
    Select and return the racemic space groups in the given sg_list
    sg_list should be a list of integers
    list of chiral space groups: 
    http://www-chimie.u-strasbg.fr/csd.doc/ConQuest/PortableHTML/conquest_portable-3-325.html
    '''
    result = []
    for sg in sg_list:
        if not sg in [1]+range(3,6)+range(16,25)+range(75,81)+range(89,99)\
            +range(143,147)+range(149,156)+range(168,174)+range(177,183)\
            +range(195,200)+range(207,215):
            result.append(sg)
    return result


class symop():
    def __init__(self,st):
        self.mat=[]
        for i in range (0,3):
            if st[i*2]=='+':
                nn=1
            else:
                nn=-1
            if st[i*2+1]=='x':
                np=0
            if st[i*2+1]=='y':
                np=1
            if st[i*2+1]=='z':
                np=2
            nline=[0,0,0]
            if st[i*2+1]!='0':
                nline[np]+=nn
            self.mat.append(nline)

