from abaqus import *
from abaqusConstants import *
import random
import math
from datetime import datetime

# Create a new model
myModel = mdb.Model(name='Model-1')

# Inclusion parameters
inclusion_axis_1 = 2.5e-3  # Major axis
inclusion_axis_2 = 2.5e-3  # Minor axis
inclusion_axis_3 = 2.5e-3  # Third axis (for ellipsoid)
inclusion_volume = 4/3 * math.pi * inclusion_axis_1 * inclusion_axis_2 * inclusion_axis_3  # Volume of ellipsoid

# Cylinder parameters
cylinder_radius = 0.05     
cylinder_height = 0.05
cylinder_volume = math.pi * cylinder_radius**2 * cylinder_height

# Set volume fraction and calculate number of inclusions
volume_fraction = 0.285
num_inclusions = int(volume_fraction * cylinder_volume / inclusion_volume)

# Create inclusion geometry
def create_inclusion_sketch():
    s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s.FixedConstraint(entity=g[2])
    s.EllipseByCenterPerimeter(center=(0.0, 0.0), 
                              axisPoint1=(0.0, inclusion_axis_1),
                              axisPoint2=(-inclusion_axis_2, inclusion_axis_1*0.75))
    s.CoincidentConstraint(entity1=v[0], entity2=g[2], addUndoState=False)
    s.CoincidentConstraint(entity1=v[2], entity2=g[2], addUndoState=False)
    s.autoTrimCurve(curve1=g[3], point1=(-2.60297012329102, -0.512653350830078))
    s.Line(point1=(0.0, inclusion_axis_1), point2=(0.0, -inclusion_axis_1))
    s.VerticalConstraint(entity=g[6], addUndoState=False)
    s.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)
    return s

# Create inclusion part
s = create_inclusion_sketch()
p = myModel.Part(name='Inclusion', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p.BaseSolidRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
s.unsetPrimaryObject()

# Create cylinder part
def create_cylinder(myModel, diameter, height):
    s = myModel.ConstrainedSketch(name='cylinder_sketch', sheetSize=200.0)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(diameter/2, 0.0))
    p = myModel.Part(name='cylinder_part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=height)
    return p

cylinder_part = create_cylinder(myModel, diameter=2*cylinder_radius, height=cylinder_height)

# Create assembly
myAssembly = myModel.rootAssembly
cylinder_instance = myAssembly.Instance(name='cylinder_instance', 
                                      part=cylinder_part, 
                                      dependent=OFF)

# Adjusted buffer factor and max attempts
buffer_factor = 1.005  # Slightly reduced buffer
max_attempts = num_inclusions * 500  # Further increased attempts

# Improved intersection check function
def does_intersect(x, y, z, inclusion_coordinates, min_distance):
    min_distance_with_buffer = min_distance * buffer_factor
    for coords in inclusion_coordinates:
        distance = math.sqrt((coords[0] - x)**2 + 
                             (coords[1] - y)**2 + 
                             (coords[2] - z)**2)
        if distance < min_distance_with_buffer:
            return True
    return False

# Random position generator
def generate_random_position(cylinder_radius, cylinder_height):
    while True:
        x = random.uniform(-cylinder_radius, cylinder_radius)
        y = random.uniform(-cylinder_radius, cylinder_radius)
        if x*x + y*y <= cylinder_radius*cylinder_radius:
            break
    z = random.uniform(-cylinder_height/2, cylinder_height/2)  # Centered around origin
    return x, y, z

# Store successful inclusion coordinates
inclusion_coordinates = []

# Minimum distance between inclusions
min_distance = 2.0 * max(inclusion_axis_1, inclusion_axis_2, inclusion_axis_3)

# Initialize counter for instance naming
instance_counter = 1

# Place inclusions with dynamic adjustment
successful_inclusions = 0
attempt_count = 0
batch_size = 10

print("Starting inclusion placement...")
while successful_inclusions < num_inclusions and attempt_count < max_attempts:
    attempt_count += 1
    
    x, y, z = generate_random_position(cylinder_radius, cylinder_height)
    
    if math.sqrt(x**2 + y**2) <= (cylinder_radius - inclusion_axis_2 * 1.1):
        if not does_intersect(x, y, z, inclusion_coordinates, min_distance):
            try:
                inclusion_part = myModel.parts['Inclusion']
                instance_name = 'Inc_Instance_%d_%s' % (instance_counter, datetime.now().strftime('%Y%m%d%H%M%S'))
                
                inclusion_instance = myAssembly.Instance(
                    name=instance_name, 
                    part=inclusion_part, 
                    dependent=OFF
                )
                
                inclusion_instance.translate(vector=(x, y, z))
                
                inclusion_coordinates.append((x, y, z))
                successful_inclusions += 1
                instance_counter += 1
                
                if successful_inclusions % batch_size == 0:
                    print("Placed %d/%d inclusions" % (successful_inclusions, num_inclusions))
                    
            except Exception as e:
                print("Failed to create inclusion %d at (%.4f, %.4f, %.4f): %s" % (instance_counter, x, y, z, str(e)))
                continue

# Create sets for each inclusion
for i in range(1, instance_counter):
    instance_name = 'Inc_Instance_%d_%s' % (i, datetime.now().strftime('%Y%m%d%H%M%S'))
    inclusion_instance = myAssembly.instances[instance_name]
    myAssembly.Set(name='Set_%s' % instance_name, instances=(inclusion_instance,))

# Create a set for the cylinder
myAssembly.Set(name='Set_cylinder', instances=(cylinder_instance,))

# Merge inclusions with the cylinder
all_instances = [myAssembly.instances['Inc_Instance_%d_%s' % (i, datetime.now().strftime('%Y%m%d%H%M%S'))] for i in range(1, instance_counter)]
all_instances.append(cylinder_instance)

merged_part = myAssembly.InstanceFromBooleanMerge(
    name='Merged_Part',
    instances=all_instances,
    keepIntersections=ON,
    originalInstances=DELETE,
    domain=GEOMETRY
)

# Create materials
def create_material(myModel, materialName, youngsModulus, poissonRatio):
    material = myModel.Material(name=materialName)
    material.Elastic(table=((youngsModulus, poissonRatio), ))

# Create matrix and inclusion materials
create_material(myModel, 'MatrixMaterial', 200.0E9, 0.29)
create_material(myModel, 'InclusionMaterial', 100.0E9, 0.3)

# Assign materials to parts
matrix_material = myModel.materials['MatrixMaterial']
inclusion_material = myModel.materials['InclusionMaterial']
cylinder_part.setMaterial(matrix_material)
inclusion_part.setMaterial(inclusion_material)

# Final statistics
print("\nPlacement Complete!")
print("Number of successful inclusions: %d" % successful_inclusions)
print("Volume of one inclusion: %.6e cubic meters" % inclusion_volume)
total_inclusion_volume = successful_inclusions * inclusion_volume
print("Total volume of inclusions: %.6e cubic meters" % total_inclusion_volume)
achieved_volume_fraction = (total_inclusion_volume / cylinder_volume) * 100
print("Achieved volume fraction: %.2f%%" % achieved_volume_fraction)
print("Target volume fraction: %.2f%%" % (volume_fraction * 100))
print("Placement attempts: %d" % attempt_count)
print("Success rate: %.2f%%" % (successful_inclusions/float(attempt_count)*100))