import numpy as np

def saver(obj,location,string_base):
    x1_size=obj[:,0,0,0].size
    x2_size=obj[0,:,0,0].size
    i_size=obj[0,0,:,0].size
    j_size=obj[0,0,0,:].size
    meta_data=np.array([x1_size,x2_size,i_size,j_size])
    
    np.savetxt(location+'meta_data_'+string_base+'.txt',meta_data.astype(int))
    for i in range(i_size):
        #print(str(i+1)+'/'+str(i_size))
        for j in range(j_size):
            string=location+string_base+'_i'+str(i)+'_j'+str(j)+'.txt'
            np.savetxt(string,obj[:,:,i,j])
    print('saved: ',string)

def loader(location,string_base):
    meta_data=np.loadtxt(location+'meta_data_'+string_base+'.txt')
    x1=int(meta_data[0])
    x2=int(meta_data[1])
    i_size=int(meta_data[2])
    j_size=int(meta_data[3])
    
    out=np.zeros((x1,x2,i_size,j_size))+0j
    
    for i in range(i_size):
        for j in range(j_size):
            string=location+string_base+'_i'+str(i)+'_j'+str(j)+'.txt'
            out[:,:,i,j]=np.loadtxt(string,dtype=complex)
    return out

def saver_3d(obj,location,string_base):
    x1_size=obj[:,0,0].size
    x2_size=obj[0,:,0].size
    i_size=obj[0,0,:].size
    meta_data=np.array([x1_size,x2_size,i_size])
    
    np.savetxt(location+'meta_data_'+string_base+'.txt',meta_data.astype(int))
    for i in range(i_size):
        #print(str(i+1)+'/'+str(i_size))
        string=location+string_base+'_i'+str(i)+'.txt'
        np.savetxt(string,obj[:,:,i])
    print('saved: ',string)


def loader_3d(location,string_base):
    meta_data=np.loadtxt(location+'meta_data_'+string_base+'.txt')
    x1=int(meta_data[0])
    x2=int(meta_data[1])
    i_size=int(meta_data[2])
    
    out=np.zeros((x1,x2,i_size))+0j
    
    for i in range(i_size):
        string=location+string_base+'_i'+str(i)+'.txt'
        out[:,:,i]=np.loadtxt(string,dtype=complex)
    return out
