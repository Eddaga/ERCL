import torch
import torch.nn as nn
import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

import matplotlib.pyplot as plt

from sklearn.datasets import load_wine
from sklearn.metrics import accuracy_score

# Define hyperparameters
rowtoskip = 2500
hidden_size = 15
output_size = 1
num_epochs = 2500
batch_size = 7000
learning_rate = 0.1
outlierClear = 0

makeModel = 1
saveModel = 0

PATH = "hello~"

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
if device == 'cpu':
    exit()
else :
    print("cuda working\r\n")
# Data NameTag
#       A
#       tick
#       B       C       D   
#       break   accel   pittot
#       E       F
#       bTemp   bPress
#       G       H       I       J       K       L
#       aAngle  aC      aTemp   aTq     aVel    aV
#       M       N       O       P       Q       R
#       bAngle  bC      bTemp   bTq     bVel    bV
#       S       T       U       V
#       fTemp   fC      fSOC    fV
#       W       X       Y       Z
#       bTemp   bC      bSOC    bV
#       AA      AB      AC      AD
#       tTemp   tC      tSOC    tV
#       AE      AF
#       latti  longi
#       AG      AH      AI      AJ      AK      AL      AM
#       AccX    AccY    AccZ    GyroX   GyroY   GyroZ   Temp
#       AN      AO      AP
#       Sec     Min     Hour

# Load data
#data = pd.read_excel('data.xlsx', header=None, usecols="B:S,U,W,Y,AA,AC,AE:AM", skiprows=rowtoskip)
#input_size = 32
#print("B:S,U,W,Y,AA,AC,AE:AM")

data = pd.read_excel('data.xlsx', header=None, usecols="C,K,AG:AL", skiprows=rowtoskip)
input_size = 8
print("C,K,AG:AL")

#data = pd.read_excel('data.xlsx', header=None, usecols="C, H, N,L,R", skiprows=rowtoskip)
#input_size = 5
#print("C H N L R")

#data = pd.read_excel('data.xlsx', header=None, usecols="B,C, H,K,L, N,Q,R ", skiprows=rowtoskip)
#input_size = 8
#print("B,C, H,K,L, N,Q,R")

#data = pd.read_excel('data.xlsx', header=None, usecols="B,C, H,K,L, N,Q,R, AG:AL", skiprows=rowtoskip)
#input_size = 14
#print("input = B,C, H,K,L, N,Q,R, AG:AL\n")

#data = pd.read_excel('data.xlsx', header=None, usecols="H,J,L,N,P,R", skiprows=rowtoskip)
#input_size = 6



target = pd.read_excel('data.xlsx', header=None, usecols="T,V,Z,AD", skiprows=rowtoskip)
target = (target.iloc[:, 3] + target.iloc[:, 2] + target.iloc[:, 1]) * target.iloc[:, 0]
print("target = (V+Z+AD) * T\n")
"""
wine = load_wine()
data = wine.data
target = wine.target
print(data)
print(target)
"""

def remove_outliers(data, target, z=3):
    """
    Remove outliers from data and target using Z-score method.
    """
    # Compute Z-score for each column in data and target
    data_zscore = np.abs(data - data.mean()) / data.std()
    target_zscore = np.abs(target - target.mean()) / target.std()
    
    # Identify rows with Z-score greater than z
    rows_to_remove = (data_zscore > z).any(axis=1) | (target_zscore > z)
    
    # Remove rows with Z-score greater than z
    data_new = data[~rows_to_remove]
    target_new = target[~rows_to_remove]
    
    return data_new, target_new

if outlierClear:
    print("data outlier Start")
    data, target = remove_outliers(data, target)


# Split data into train, validation, test sets
train_data, test_data, train_target, test_target = train_test_split(data, target, test_size=0.2, random_state=42)
train_data, val_data, train_target, val_target = train_test_split(train_data, train_target, test_size=0.25, random_state=42)
print(len(train_data))
# Normalize data
input_scaler = MinMaxScaler()
train_data = input_scaler.fit_transform(train_data)
test_data = input_scaler.transform(test_data)
val_data = input_scaler.transform(val_data)

train_target = train_target.values.reshape(-1,1)
test_target = test_target.values.reshape(-1,1)
val_target = val_target.values.reshape(-1,1)

# Define the neural network
class MLP(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, output_size)
    
    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        return out


# Move data to device
train_data = torch.tensor(train_data, dtype=torch.float32).to(device)
train_target = torch.tensor(train_target, dtype=torch.float32).to(device)
val_data = torch.tensor(val_data, dtype=torch.float32).to(device)
val_target = torch.tensor(val_target, dtype=torch.float32).to(device)
test_data = torch.tensor(test_data, dtype=torch.float32).to(device)
test_target = torch.tensor(test_target, dtype=torch.float32).to(device)


# Train the model

if makeModel:
    model = MLP(input_size, hidden_size, output_size).to(device)
else:
    model = torch.load(PATH).to(device)
    model.eval()

criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

total_step = len(train_data) // batch_size
for epoch in range(num_epochs):
    for i in range(total_step):
        # Obtain a batch of training data
        batch_data = train_data[i*batch_size:(i+1)*batch_size]
        batch_target = train_target[i*batch_size:(i+1)*batch_size]

        # Forward pass
        outputs = model(batch_data)
        loss = criterion(outputs, batch_target)

        # Backward and optimize
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    # Print training loss for each epoch
    print('Epoch [{}/{}], Loss: {:.4f}'.format(epoch+1, num_epochs, loss.item()))


# Evaluate the model
model.eval()
with torch.no_grad():
        # Compute validation loss
    val_outputs = model(val_data)
    val_loss = criterion(val_outputs, val_target)
    print('Validation Loss: {:.4f}'.format(val_loss.item()))
    

    # Compute test loss
    test_outputs = model(test_data)
    test_loss = criterion(test_outputs, test_target)
    print('Test Loss: {:.4f}'.format(test_loss.item()))
    
if saveModel:
    torch.save(model.state_dict(), PATH)

def write_dataframes_to_excel(train_target, test_target, val_target, outputs, val_outputs, test_outputs):
    with pd.ExcelWriter('dataframes.xlsx') as writer:
        train_target.to_excel(writer, sheet_name='train_target')
        test_target.to_excel(writer, sheet_name='test_target')
        val_target.to_excel(writer, sheet_name='val_target')
        outputs.to_excel(writer, sheet_name='train_outputs')
        test_outputs.to_excel(writer, sheet_name='test_outputs')
        val_outputs.to_excel(writer, sheet_name='val_outputs')

print("xlsx writng start\r\n")
train_target = pd.DataFrame(train_target.detach().cpu().numpy())
test_target = pd.DataFrame(test_target.detach().cpu().numpy())
val_target = pd.DataFrame(val_target.detach().cpu().numpy())
outputs = pd.DataFrame(outputs.detach().cpu().numpy())
test_outputs = pd.DataFrame(test_outputs.detach().cpu().numpy())
val_outputs = pd.DataFrame(val_outputs.detach().cpu().numpy())

write_dataframes_to_excel(train_target, test_target, val_target, outputs, val_outputs, test_outputs)
print("xlsx writng done\r\n")