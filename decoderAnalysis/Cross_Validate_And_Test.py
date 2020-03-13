import tensorflow as tf
import keras
import numpy as np
from keras.layers import Input, Dense, Dropout
from keras.layers import BatchNormalization
from keras.models import Model
from keras import regularizers
import keras.backend as K
from keras.losses import mean_squared_error
import pickle
import keras
import json

from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
# RMSE loss definition
def RMSE_loss(y_pred,y_true):
    loss=np.sum((y_pred-y_true)**2,axis=-1)
    loss=loss/loss.shape[0]
    loss=np.sum(loss)
    return(np.sqrt(loss))
  
# Code for encoding data as JSON    
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
# Function for loading and testing cross validation scores on decoded data
def Cross_Validation_Testing(base_neuron_path,base_neuron_position,test_neuron_paths=None,test_position_paths=None,neuron_subset=None):
    #Input parameters for cross validoation testing
    Result_Data={};
    Result_Data['Parameters']={};
    Result_Data['Parameters']['num_folds']=10 #number of cross validation folds
    Result_Data['Parameters']['window_length']=20 #Window length for determining how many time instance of data are used as input to the network
    Result_Data['Parameters']['offsets']=4 # time offset in frames between behavior recording and neural recording. Previous research indicated a latency of approximately 250 ms
    Result_Data['Parameters']['num_neurons']=700 #Number of neurons included in the first layer of the neural network
    Result_Data['Parameters']['num_layers']=5 #Number of layers in the network
    
    Result_Data['Parameters']['iterations']=700 #Number of epochs to run data. 
    Result_Data['Parameters']['cutoff']=.1 #For linear track data, this removes data corresponding to time points where the mouse is at the end of the track. If this value is .1, 10% is removed from both ends
    Result_Data['Parameters']['tracklength']=100 #Length of the linear track (in centimeters). Track lenght is normalized to this value as our data was frequently scaled incorrectly initially
    
    
    Result_Data['neural_Data']=[]
    #Neural data input loading. Assumes comma delimited or default numpy delimited data
    try:
        Result_Data['neural_Data'].append(np.transpose(np.loadtxt(base_neuron_path,delimiter=','))) 
    except: 
        Result_Data['neural_Data'].append(np.transpose(np.loadtxt(base_neuron_path)))
    Result_Data['position_Data']=[]
    #Loading position data
    try:
        Result_Data['position_Data'].append(np.loadtxt(base_neuron_position,delimiter=','))
    except:
        Result_Data['position_Data'].append(np.loadtxt(base_neuron_position))
    #Construd standard scaler (see scikit learn) and apply to neural data
    stand=StandardScaler(with_mean=False)
    Result_Data['neural_Data'][0]=stand.fit_transform(Result_Data['neural_Data'][0])



    #Thresholds neural data at 10% of the maximum, and zeros out all data with values below this threshold (for noise reduction)
    maxim=.1*np.max(Result_Data['neural_Data'][0],axis=0)
    for i in range(Result_Data['neural_Data'][0].shape[1]):
        indices=Result_Data['neural_Data'][0][:,i]<maxim[i]
        Result_Data['neural_Data'][0][indices,i]=0
    for i in range(0,len(test_position_paths)):
    #Loads any test data based on paths submitted to the function. Repeats preprocessing
        try:
            Result_Data['neural_Data'].append(np.transpose(np.loadtxt(test_neuron_paths[i],delimiter=',')))
        except:
            Result_Data['neural_Data'].append(np.transpose(np.loadtxt(test_neuron_paths[i])))
        stand=StandardScaler(with_mean=False)    
        Result_Data['neural_Data'][-1]=stand.fit_transform(Result_Data['neural_Data'][-1])
        
        maxim=.1*np.max(Result_Data['neural_Data'][-1],axis=0)
        for j in range(Result_Data['neural_Data'][-1].shape[1]):
            indices=Result_Data['neural_Data'][-1][:,j]<maxim[j]
            Result_Data['neural_Data'][-1][indices,i]=0
        
        try:
            Result_Data['position_Data'].append(np.loadtxt(test_position_paths[i]))
        except:
            Result_Data['position_Data'].append(np.loadtxt(test_position_paths[i],delimiter=','))
    
    Result_Data['Parameters']['scale']=[];
    #Determine scaling value for position data so that tracklength is correct
    for i in range(len(Result_Data['position_Data'])):
        
        Result_Data['Parameters']['scale'].append(Result_Data['Parameters']['tracklength']/np.max(Result_Data['position_Data'][i]))
        
    
    
    
    
    #Collects neuron subset into Results, if a subset of neurons is used  (such as indices of place cells) This is provided as the path to a txt file
    Result_Data['Parameters']['neuron_subset']=neuron_subset
   
    if Result_Data['Parameters']['neuron_subset'] is not None:
        neuron_subset=Result_Data['Parameters']['neuron_subset']
        try:
            neuron_subset=np.loadtxt(neuron_subset)
            #Currently this next line removes 1 from each index, as indices were calculated using matlab, which counts from 1 instead of 0
            neuron_subset=neuron_subset.astype('int')-1
            Result_Data['Parameters']['neuron_subset']=neuron_subset
        except:
            pass
    	#Restrict to neurons in the subset
        for i in range(len(Result_Data['neural_Data'])): 
            Result_Data['neural_Data'][i]=Result_Data['neural_Data'][i][:,neuron_subset]
    
    
    
    
    for i in range(len(Result_Data['neural_Data'])):
       	#Use specified window to construct neural network input of size (number of data points)x((2*window_length+1)*number of neurons)
        reshaped_Data=np.zeros((Result_Data['neural_Data'][i].shape[0]-2*Result_Data['Parameters']['window_length'],(2*Result_Data['Parameters']['window_length']+1)*Result_Data['neural_Data'][i].shape[1]))
        new_positions=np.zeros((reshaped_Data.shape[0],2))
        for j in range(Result_Data['neural_Data'][i].shape[0]-2*Result_Data['Parameters']['window_length']):
            reshaped_Data[j,:]=np.reshape(Result_Data['neural_Data'][i][j:j+2*Result_Data['Parameters']['window_length']+1,:],[1,-1])
            new_positions[j,:]=Result_Data['position_Data'][i][j+Result_Data['Parameters']['window_length']+1,:]
            
        new_positions=new_positions*Result_Data['Parameters']['scale'][i]
        
        #new_positions=new_positions[220:,:]
        #reshaped_Data=reshaped_Data[220:,:]
        offset=Result_Data['Parameters']['offsets']
        reshaped_Data=reshaped_Data[offset+1:,:]
        new_positions=new_positions[:-offset-1:]
       
        #Remove Data based on cutoff
        if Result_Data['Parameters']['cutoff']!=0:
            indices=(new_positions[:,0]>Result_Data['Parameters']['cutoff']* Result_Data['Parameters']['tracklength'])&(new_positions[:,0]<Result_Data['Parameters']['tracklength']-Result_Data['Parameters']['cutoff']* Result_Data['Parameters']['tracklength'])
            new_positions=new_positions[indices,:]
            #new_positions[:,0]=new_positions[:,0]-Result_Data['Parameters']['tracklength']/2
            reshaped_Data=reshaped_Data[indices,:]
        
        Result_Data['neural_Data'][i]=reshaped_Data
        new_positions[:,0]=new_positions[:,0]-Result_Data['Parameters']['tracklength']/2
        
        Result_Data['position_Data'][i]=new_positions
    #Initialize results as empty sets    
    Result_Data['train_neural'] =[]
    Result_Data['test_neural']=[]
    Result_Data['train_position']=[]
    Result_Data['test_position']=[]
    Result_Data['train_predictions']=[]
    Result_Data['test_predictions']=[]
    Result_Data['secondary_predictions']=[[] for i in range(len(test_neuron_paths))]
    
    Result_Data['train_R2']=[]
    Result_Data['test_R2']=[]
    Result_Data['secondary_R2']=[[] for i in range(len(test_neuron_paths))]
    
    Result_Data['train_RMSE']=[]
    Result_Data['test_RMSE']=[]
    Result_Data['secondary_RMSE']=[[] for i in range(len(test_neuron_paths))]
    #Create train and test sets and normalize using standard scaler
    for j in range(Result_Data['Parameters']['num_folds']):
        test_indices=np.arange((int(j*Result_Data['neural_Data'][0].shape[0]/Result_Data['Parameters']['num_folds'])),
                                   int((j+1)*Result_Data['neural_Data'][0].shape[0]/Result_Data['Parameters']['num_folds']))
        train_indices=np.setdiff1d(np.arange(Result_Data['neural_Data'][0].shape[0]),test_indices)
        
        train_neural=Result_Data['neural_Data'][0][train_indices,:]
        test_neural=Result_Data['neural_Data'][0][test_indices,:]
        stand=StandardScaler()
        train_neural=stand.fit_transform(train_neural)
        test_neural=stand.transform(test_neural)
        
        
       
        train_position=Result_Data['position_Data'][0][train_indices,:]
        test_position=Result_Data['position_Data'][0][test_indices,:]
        
        
        
        Result_Data['train_neural'].append(train_neural)
        Result_Data['test_neural'].append(test_neural)
        Result_Data['train_position'].append(train_position)
        Result_Data['test_position'].append(test_position)
        
        K.clear_session()
        # Construct neural network   
        input_matrix=Input(shape=(train_neural.shape[1],))
        encoded=Dense(train_neural.shape[1],activation='relu',kernel_initializer=keras.initializers.glorot_normal(seed=None),kernel_regularizer=regularizers.l2(.05))(input_matrix)
        batch=BatchNormalization()(encoded)
        drop=Dropout(.1)(batch)
        for k in range(Result_Data['Parameters']['num_layers']):
            if int(Result_Data['Parameters']['num_neurons']/((2)**k))>2:
		#Number of neurons per layer decreases by a factor of 2 per layer
                encoded=Dense(int(Result_Data['Parameters']['num_neurons']/((2)**k)),activation='relu',kernel_initializer=keras.initializers.glorot_normal(seed=None),kernel_regularizer=regularizers.l2(.05))(drop)
                
                batch=BatchNormalization()(encoded)
                drop=Dropout(.1)(batch)
                
            else:
                break
            
        decoded=Dense(2)(drop)
        position_decoder=Model(input_matrix,decoded)
        #Train network
        position_decoder.compile(optimizer='adam',loss='mean_squared_error')
        history=position_decoder.fit(train_neural,train_position,batch_size=1000,callbacks=[keras.callbacks.EarlyStopping(patience=100,restore_best_weights=True)],epochs=Result_Data['Parameters']['iterations'],validation_split=0,
                                    validation_data=(test_neural,test_position),shuffle=True)
	#Record predictions, calculat R2 and RMSE metrics        
	Result_Data['train_predictions'].append(position_decoder.predict(train_neural))
        Result_Data['test_predictions'].append(position_decoder.predict(test_neural))
            
        Result_Data['train_R2'].append(r2_score(Result_Data['train_predictions'][-1],Result_Data['train_position'][-1]))
        Result_Data['test_R2'].append(r2_score(Result_Data['test_predictions'][-1],Result_Data['test_position'][-1]))
            
        Result_Data['train_RMSE'].append(RMSE_loss(Result_Data['train_predictions'][-1],Result_Data['train_position'][-1]))
        Result_Data['test_RMSE'].append(RMSE_loss(Result_Data['test_predictions'][-1],Result_Data['test_position'][-1]))
            
            
            
            
        for l in range(len(test_neuron_paths)):
            Result_Data['secondary_predictions'][l].append(position_decoder.predict(Result_Data['neural_Data'][l+1]))
            Result_Data['secondary_R2'][l].append(r2_score( Result_Data['secondary_predictions'][l][-1],Result_Data['position_Data'][l+1]))
            Result_Data['secondary_RMSE'][l].append(RMSE_loss( Result_Data['secondary_predictions'][l][-1],Result_Data['position_Data'][l+1]))
        K.clear_session()
        del position_decoder
        tf.reset_default_graph()
    #Save result data in appropriate location
    with open('M3244_place.json','w') as fp:
        json.dump(Result_Data,fp,cls=NumpyEncoder)
"""
secondary=['/home/user1/Desktop/Research/LinearTrackData/M3243_new/CNO_neuron.txt',
           '/home/user1/Desktop/Research/LinearTrackData/M3243_new/PCTRL_neuron.txt']
test_position_paths=['/home/user1/Desktop/Research/LinearTrackData/M3243_new/CNO_pos.txt',
                     '/home/user1/Desktop/Research/LinearTrackData/M3243_new/PCTRL_pos.txt'  ]
base_neuron_path='/home/user1/Desktop/Research/LinearTrackData/M3243_new/CTRL_neuron.txt'
base_neuron_position='/home/user1/Desktop/Research/LinearTrackData/M3243_new/CTRL_pos.txt'  
neuron_subset='/home/user1/Desktop/Research/LinearTrackData/M3243_new/PlaceCells.txt'      


"""
#Perform cross validation on data as shown. All data was saved as .txt files prior to analysis.
base_neuron_path='/home/user1/Desktop/Research/LinearTrackData/M3244/0217/0217C_noise_scaled.txt'
base_neuron_position='/home/user1/Desktop/Research/LinearTrackData/M3244/0217/0217pos.txt'

secondary=[
           '/home/user1/Desktop/Research/LinearTrackData/M3244/0221/0217C_noise_scaled.txt',
          '/home/user1/Desktop/Research/LinearTrackData/M3244/0224/0217C_noise_scaled.txt'
          
          ]

test_position_paths=['/home/user1/Desktop/Research/LinearTrackData/M3244/0221/0217pos.txt',
          '/home/user1/Desktop/Research/LinearTrackData/M3244/0224/0217pos.txt'
         ]
neuron_subset='/home/user1/Desktop/Research/LinearTrackData/M3244/commonplacecells.txt'

#neuron_subset=None (Delete this if no subset is desired)


#Run cross validation
Cross_Validation_Testing(base_neuron_path,base_neuron_position,secondary,test_position_paths,neuron_subset=neuron_subset)
            
            

            
            
            
    
        
        
        
        
        
    
    
    
    
    
    
    
    
    



    




