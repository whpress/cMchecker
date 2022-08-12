TEST(BadValuesTest, BadFatherSonRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/bad/father_son.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE((returned_chisquared > 0) && (returned_chisquared < 0.1));
}

TEST(BadValuesTest, BadUncleNephewRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/bad/uncle_nephew.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE((returned_chisquared > 0) && (returned_chisquared < 0.1));
}

TEST(BadValuesTest, BadSiblingRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/bad/sibling.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE((returned_chisquared > 0) && (returned_chisquared < 0.1));
}

TEST(BadValuesTest, BadHalfSiblingRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/bad/half_sibling.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE((returned_chisquared > 0) && (returned_chisquared < 0.1));
}
