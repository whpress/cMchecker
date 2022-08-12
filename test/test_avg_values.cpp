TEST(AvgValuesTest, AvgFatherSonRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/avg/father_son.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE((returned_chisquared > 0.6) && (returned_chisquared < 0.7));
}

TEST(AvgValuesTest, AvgUncleNephewRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/avg/uncle_nephew.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE((returned_chisquared > 0.8) && (returned_chisquared < 0.9));
}

TEST(AvgValuesTest, AvgSiblingRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/avg/sibling.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE((returned_chisquared > 0.009) && (returned_chisquared < 0.02));
}

TEST(AvgValuesTest, AvgHalfSiblingRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/avg/half_sibling.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE((returned_chisquared > 0.75) && (returned_chisquared < 0.85));
}
