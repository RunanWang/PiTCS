#include <filesystem>
#include <fmt/format.h>
#include "gflags/gflags.h"
#include "glog/logging.h"

#include "Utils/Timer.h"
#include "Utils/Memory.h"
#include "Utils/logging.h"
#include "DataStructure/MLGWithSchema.h"
#include "test/testHeader.h"
#include "Algorithm/PivotalTruss.h"
#include "constant.h"

using namespace std;

DEFINE_string(output_path, "./output", "output dir");
DEFINE_string(input_path, "./Dataset", "input graph file");
DEFINE_string(dataset, "RM", "input dataset name");
DEFINE_string(algo, "PiTCS-Index", "execute algorithm. PiTCS-Index/RPiTCS-Index/PiTCS-online");
DEFINE_string(datatype, "label", "Dataset type. label/no-label");

vector<NODE_TYPE> queries = {
        310, 326, 659, 676, 772,
};

void check_dir(const filesystem::path &path) {
    if (not filesystem::exists(path)) {
        filesystem::create_directory(path);
    }
}

void routinePivot() {
    auto output_path = filesystem::path(FLAGS_output_path);
    auto input_file_path = filesystem::path(FLAGS_input_path).append(fmt::format("{}.txt", FLAGS_dataset));
    auto input_path = filesystem::path(FLAGS_input_path);
    check_dir(output_path);
    // load graph


    if (FLAGS_algo == "decomposition") {
        LOG(INFO) << fmt::format("Reading Graph from {}...", input_file_path.string());
        auto timer = new Timer();
        timer->startTimer();
        auto mlg = new MLGWithSchema();
        if (not filesystem::exists(input_file_path)) {
            LOG(ERROR) << fmt::format("No Such File: {}", input_file_path.string());
            return;
        }
        mlg->LoadFromFile(input_file_path);
        LOG(INFO) << fmt::format("Graph Size = {:.4f} MB", mlg->getMemUsage());
        timer->endTimer();
        LOG(INFO) << fmt::format("Loading Phase: {:.4f} Seconds.", timer->getTimerSecond());
        LOG(INFO) << fmt::format("Peak Memory Usage in loading = {:.4f} MB", GetPeakRSSInMB());
        LAYER_TYPE maxComboSize = ALLOW_NO_PIVOT
                                  ? mlg->layersNum * (mlg->layersNum - 1) + mlg->layersNum
                                  : mlg->layersNum * (mlg->layersNum - 1);
        auto result = new EDGE_TYPE *[maxComboSize];
        auto ci = new PivotalComboIndex(mlg->layersNum);
        PivotTrussDecomposition(mlg, result, ci, ALLOW_NO_PIVOT);
        for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
            delete[] result[i];
        }
        delete[] result;
        delete ci;
    } else if (FLAGS_algo == "index") {
        LOG(INFO) << fmt::format("Reading Graph from {}...", input_file_path.string());
        auto timer = new Timer();
        timer->startTimer();
        auto mlg = new MLGWithSchema();
        if (not filesystem::exists(input_file_path)) {
            LOG(ERROR) << fmt::format("No Such File: {}", input_file_path.string());
            return;
        }
        mlg->LoadFromFile(input_file_path);
        LOG(INFO) << fmt::format("Graph Size = {:.4f} MB", mlg->getMemUsage());
        timer->endTimer();
        LOG(INFO) << fmt::format("Loading Phase: {:.4f} Seconds.", timer->getTimerSecond());
        LOG(INFO) << fmt::format("Peak Memory Usage in loading = {:.4f} MB", GetPeakRSSInMB());
        LAYER_TYPE maxComboSize = ALLOW_NO_PIVOT
                                  ? mlg->layersNum * (mlg->layersNum - 1) + mlg->layersNum
                                  : mlg->layersNum * (mlg->layersNum - 1);
        LOG(INFO) << "### Constructing PivotTruss Index ###";
        auto trussness = new EDGE_TYPE *[maxComboSize];
        auto ci = new PivotalComboIndex(mlg->layersNum);
        PivotTrussDecomposition(mlg, trussness, ci, ALLOW_NO_PIVOT);
        auto **index = new EqualTree *[ci->comboSize + 1];
        PivotIndexAll(mlg, index, ci, trussness);
        double totalIndexSize = 0;
        for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
            totalIndexSize += index[i]->getMemUsage();
            delete index[i];
        }
        delete ci;
    } else if (FLAGS_algo == "PiTCS-Index") {
        if (FLAGS_datatype == "label") {
            if (FLAGS_dataset == "RM") {
                LOG(INFO) << fmt::format("[Exp=SOTA-GT-RM] Generating density and f1 of community for RM.");
                testPiRMCommunity(input_path);
            } else if (FLAGS_dataset == "terrorist") {
                LOG(INFO)
                        << fmt::format("[Exp=SOTA-GT-terrorist] Generating density and f1 of community for terrorist.");
                testPiTerroristCommunity(input_path);
            } else {
                LOG(ERROR) << fmt::format("Dataset {} dose not have ground truth community.", FLAGS_dataset);
            }
        } else {
            LOG(INFO) << fmt::format("[Exp=SOTA-GIVEN-PIVOT] Generating density and diameter of community for {}.",
                                     FLAGS_dataset);
            testPiCommunityByDensityUsingGivenQueries(input_path, FLAGS_dataset, queries);
        }
    } else if (FLAGS_algo == "RPiTCS-Index") {
        if (FLAGS_datatype == "label") {
            if (FLAGS_dataset == "RM") {
                LOG(INFO) << fmt::format("[Exp=SOTA-GT-RM] Generating density and f1 of community for RM.");
                testPiWRMCommunity(input_path);
            } else if (FLAGS_dataset == "terrorist") {
                LOG(INFO)
                        << fmt::format("[Exp=SOTA-GT-terrorist] Generating density and f1 of community for terrorist.");
                testPiWTerroristCommunity(input_path);
            } else {
                LOG(ERROR) << fmt::format("Dataset {} dose not have ground truth community.", FLAGS_dataset);
            }
        } else {
            LOG(INFO) << fmt::format("[Exp=SOTA-GIVEN-RANK] Generating density and diameter of community for {}.",
                                     FLAGS_dataset);
            testRankCommunityByDensityUsingGivenQueries(input_path, FLAGS_dataset, queries);
        }
    } else if (FLAGS_algo == "PiTCS-online") {
        LOG(INFO) << fmt::format("[Exp=SOTA-NOINDEX] Generating density of community for {}.", FLAGS_dataset);
        searchPiCommByTrussnessForGivenQueries(input_path, FLAGS_dataset, queries);
    }
    // delete timer;
    // delete mlg;
}

int main(int argc, char *argv[]) {
    // init
    FLAGS_alsologtostderr = true; //设置日志消息除了日志文件之外是否去标准输出
    FLAGS_log_prefix = true; //设置日志前缀是否应该添加到每行输出
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    FLAGS_log_dir = FLAGS_output_path;

    routinePivot();

    LOG(INFO) << fmt::format("Peak Memory Usage in total = {:.4f} MB", GetPeakRSSInMB());
    return 0;
}
