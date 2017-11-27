# short parameter list for testing purposes
experiments_test <- list(
  list(initTopicAssignments = T, iter = 1),
  list(initTopicAssignments = F, iter = 1),
  list(initTopicAssignments = "SEED", iter = 1),
  list(initTopicAssignments = T, iter = 10),
  list(initTopicAssignments = F, iter = 10),
  list(initTopicAssignments = "SEED", iter = 10)
)

# finally tested configurations
experiments_final <- list(
  list(initTopicAssignments = T, iter = 1),
  list(initTopicAssignments = F, iter = 1),
  list(initTopicAssignments = "SEED", iter = 1),
  list(initTopicAssignments = T, iter = 50),
  list(initTopicAssignments = F, iter = 50),
  list(initTopicAssignments = "SEED", iter = 50),
  list(initTopicAssignments = T, iter = 100),
  list(initTopicAssignments = F, iter = 100),
  list(initTopicAssignments = "SEED", iter = 100),
  list(initTopicAssignments = T, iter = 200),
  list(initTopicAssignments = F, iter = 200),
  list(initTopicAssignments = "SEED", iter = 200),
  list(initTopicAssignments = T, iter = 300),
  list(initTopicAssignments = F, iter = 300),
  list(initTopicAssignments = "SEED", iter = 300),
  list(initTopicAssignments = T, iter = 400),
  list(initTopicAssignments = F, iter = 400),
  list(initTopicAssignments = "SEED", iter = 400),
  list(initTopicAssignments = T, iter = 500),
  list(initTopicAssignments = F, iter = 500),
  list(initTopicAssignments = "SEED", iter = 500),
  list(initTopicAssignments = T, iter = 750),
  list(initTopicAssignments = F, iter = 750),
  list(initTopicAssignments = "SEED", iter = 750),
  list(initTopicAssignments = T, iter = 1000),
  list(initTopicAssignments = F, iter = 1000),
  list(initTopicAssignments = "SEED", iter = 1000),
  list(initTopicAssignments = T, iter = 1500),
  list(initTopicAssignments = F, iter = 1500),
  list(initTopicAssignments = "SEED", iter = 1500),
  list(initTopicAssignments = T, iter = 2000),
  list(initTopicAssignments = F, iter = 2000),
  list(initTopicAssignments = "SEED", iter = 2000)
)

# set this to experiments_final to run all iteration steps
# this will take long time!
experiments <- experiments_test