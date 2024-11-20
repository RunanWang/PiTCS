//
// Created by 王润安 on 2023/12/11.
//

#ifndef MLGWORK_TESTHEADER_H
#define MLGWORK_TESTHEADER_H

#include "filesystem"
#include "../DataStructure/MLGWithSchema.h"


void testPiRMCommunity(const filesystem::path &input_dir_path);

void testPiWRMCommunity(const filesystem::path &input_dir_path);

void testPiTerroristCommunity(const filesystem::path &input_dir_path);

void testPiWTerroristCommunity(const filesystem::path &input_dir_path);

void testRankCommunityByDensityUsingGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                               const vector<NODE_TYPE> &queries);

void searchPiCommByTrussnessForGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                            const vector<NODE_TYPE> &queries);

void testPiCommunityByDensityUsingGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                               const vector<NODE_TYPE> &queries);


#endif //MLGWORK_TESTHEADER_H
