TEST(IdealValuesTest, IdealFatherSonRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/ideal/father_son.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE(returned_chisquared > 0.99);
}

TEST(IdealValuesTest, IdealUncleNephewRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/ideal/uncle_nephew.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE(returned_chisquared > 0.9);
}

TEST(IdealValuesTest, IdealSiblingRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/ideal/sibling.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE(returned_chisquared > 0.9);
}

TEST(IdealValuesTest, IdealHalfSiblingRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/ideal/half_sibling.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE(returned_chisquared > 0.9);
}