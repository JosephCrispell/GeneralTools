#!usr/usr/bin/python

# Create an array
array = [1,2,3,4,5]

# Print each value in the array
for value in array:
	print(value)

print("------------------------")

# Print each value in the array using its index
for index in range(0, len(array), 1):
	print("index: " + str(index) + " value: " + str(array[index]))

print("------------------------")	
	
# Print each value in the array in reverse using negative indices
for index in range(-1 , -len(array) - 1, -1):
	print("index: " + str(index) + " value: " + str(array[index]))