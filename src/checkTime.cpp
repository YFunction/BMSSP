#include "../include/Graph.h"
#include "../include/Dijkstra.h"
#include "../include/BMSSP.h"
#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <unordered_set>
#include <ctime>
#include <random>
#include <fstream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <atomic>
#include <vector>
#include <queue>

#ifdef _WIN64
#include <windows.h>
#endif

using namespace std;

static Graph rand_graph(int n, int outdeg, uint32_t seed){
    Graph g(n);
    // 构造固定出度 outdeg 的随机有向图；权值为 (0,1] 实数
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> vid(0, n - 1);
    std::uniform_real_distribution<double> w(0.0, 1.0);
    for (int u = 0; u < n; ++u)
    {
        std::unordered_set<int> used;
        used.insert(u);
        for (int k = 0; k < outdeg; ++k)
        {
            int v = vid(rng);
            if (v == u || used.count(v))
            {
                v = (v + 1) % n;
            }
            used.insert(v);
            g.addEdge(u, v, w(rng) + 1e-9); // break ties minimally
        }
    }
    return g;
}

/*
Graph rand_graph(int n, int deg, uint32_t seed) {
    if (n <= 0) return Graph(0);
    
    // 设置随机数生成器
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> weight_dist(0.1, 10.0); // 边权范围
    
    Graph graph(n);
    
    // 特殊情况处理
    if (deg <= 0 || n == 1) {
        return graph;
    }
    
    // 确保源点s=1（索引0）能够到达尽可能多的点
    // 使用BFS思想生成可达的图结构
    
    // 记录每个顶点已经添加的出边数量
    std::vector<int> out_degree(n, 0);
    
    // 记录每个顶点是否已被访问（从源点可达）
    std::vector<bool> visited(n, false);
    
    // 使用队列进行BFS式的图生成
    std::queue<int> q;
    q.push(0); // 源点s=1对应索引0
    visited[0] = true;
    
    // 用于记录哪些顶点已被添加到队列中
    std::vector<int> vertex_list;
    vertex_list.reserve(n);
    vertex_list.push_back(0);
    
    // 生成一个随机排列，用于随机选择目标顶点
    std::vector<int> all_vertices(n);
    for (int i = 0; i < n; ++i) all_vertices[i] = i;
    std::shuffle(all_vertices.begin(), all_vertices.end(), rng);
    
    // 优先确保图的连通性
    while (!q.empty() && vertex_list.size() < n) {
        int u = q.front();
        q.pop();
        
        // 为当前顶点添加出边，但不超过deg条
        int edges_to_add = std::min(deg - out_degree[u], 
                                   static_cast<int>(n - vertex_list.size()));
        
        if (edges_to_add > 0) {
            // 随机选择目标顶点，优先选择未访问过的
            int edges_added = 0;
            
            for (int v : all_vertices) {
                if (edges_added >= edges_to_add) break;
                
                // 跳过自己
                if (u == v) continue;
                
                // 随机决定是否添加这条边
                if (std::uniform_real_distribution<>(0, 1)(rng) < 0.7) {
                    double w = weight_dist(rng);
                    graph.addEdge(u, v, w);
                    out_degree[u]++;
                    edges_added++;
                    
                    // 如果目标顶点未访问过，将其加入队列
                    if (!visited[v]) {
                        visited[v] = true;
                        q.push(v);
                        vertex_list.push_back(v);
                    }
                }
            }
        }
        
        // 如果队列为空但还有未访问的顶点，随机选择一个加入队列
        if (q.empty() && vertex_list.size() < n) {
            for (int v : all_vertices) {
                if (!visited[v]) {
                    // 从已访问的顶点中随机选择一个作为其前驱
                    int random_predecessor = vertex_list[std::uniform_int_distribution<>(0, vertex_list.size()-1)(rng)];
                    
                    // 添加一条边
                    double w = weight_dist(rng);
                    graph.addEdge(random_predecessor, v, w);
                    out_degree[random_predecessor]++;
                    
                    visited[v] = true;
                    q.push(v);
                    vertex_list.push_back(v);
                    break;
                }
            }
        }
    }
    
    // 第二阶段：为每个顶点添加额外的随机边，使其出度接近deg
    // 但不超过deg限制
    for (int u = 0; u < n; ++u) {
        if (out_degree[u] >= deg) continue;
        
        // 还需要添加的边数
        int remaining_edges = deg - out_degree[u];
        
        // 随机选择目标顶点
        std::vector<int> candidate_vertices;
        candidate_vertices.reserve(n - 1);
        
        for (int v = 0; v < n; ++v) {
            if (u != v) candidate_vertices.push_back(v);
        }
        
        std::shuffle(candidate_vertices.begin(), candidate_vertices.end(), rng);
        
        // 添加边
        int edges_added = 0;
        for (int v : candidate_vertices) {
            if (edges_added >= remaining_edges) break;
            
            // 随机决定是否添加这条边（避免过度连接）
            if (std::uniform_real_distribution<>(0, 1)(rng) < 0.3) {
                double w = weight_dist(rng);
                graph.addEdge(u, v, w);
                out_degree[u]++;
                edges_added++;
            }
        }
    }
    
    // 第三阶段：确保图的稀疏性，如果某些顶点出度超过deg，随机删除一些边
    // 但在Graph结构中删除边比较麻烦，我们重新构建图
    Graph final_graph(n);
    
    std::uniform_real_distribution<double> final_weight_dist(0.1, 10.0);
    
    for (int u = 0; u < n; ++u) {
        // 获取原图中的所有出边
        const auto& edges = graph.adj[u];
        
        if (edges.size() <= deg) {
            // 如果不超过deg，直接复制所有边
            for (const auto& edge : edges) {
                final_graph.addEdge(u, edge.to, edge.w);
            }
        } else {
            // 如果超过deg，随机选择deg条边
            std::vector<int> indices(edges.size());
            for (size_t i = 0; i < edges.size(); ++i) indices[i] = i;
            std::shuffle(indices.begin(), indices.end(), rng);
            
            for (int i = 0; i < deg && i < static_cast<int>(edges.size()); ++i) {
                const auto& edge = edges[indices[i]];
                final_graph.addEdge(u, edge.to, edge.w);
            }
        }
    }
    
    return final_graph;
}
*/
/*
int main(){
    ios::sync_with_stdio(false);
    srand(time(NULL));
    fstream out("time.csv",ios::out|ios::app);
    out.tie(0);
    out<<"x,y\n";
    for(int round=10000;round<=10000*1000;round+=10000){
        int n = round;
        int mdeg = 3;
        int s = 0;
        uint32_t seed = rand();

        Graph g = rand_graph(n, mdeg, seed);

        auto lg = [](double x){ 
            return std::log(std::max(2.0, x)); 
        };

        double L = lg((double)n);
        int Sigma = std::max(1, (int)std::floor(std::pow(L, 1.0 / 3.0)));
        int Tau = std::max(1, (int)std::floor(std::pow(L, 2.0 / 3.0)));
        int ellTop = std::max(1, (int)std::ceil(L / std::max(1, Tau)));
        
        auto t1 = std::chrono::steady_clock::now();
        BMSSPAlgo bm(g, s, g.n, Sigma, Tau);
        BMSSPResult res;
        bool bm_ok = true;
        std::string bm_err;
        try
        {
            res = bm.run(ellTop, INF_DIST, std::vector<int>{s});
        }
        catch (const std::exception &ex)
        {
            bm_ok = false;
            bm_err = ex.what();
        }
        auto t2 = std::chrono::steady_clock::now();

        double tb = std::chrono::duration<double>(t2 - t1).count();

        out<<fixed<<setprecision(6);

        out<<n/10000<<" "<<tb*10000<<endl;

        auto d=bm.distances();

        for(auto [to,w]:g.adj[0]){
            cout<<to<<" "<<w<<endl;
            cout<<d[to]<<endl;
        }
        for(int i=0;i<d.size();i++){
            if(d[i]!=INF_DIST)cout<<i<<endl;
        }

    }
    

    return 0;
}*/

mutex g_mutex;

// 任务结构体，包含要测试的图规模
struct TestTask {
    int n;           // 顶点数
    int round;       // round值
    uint32_t seed;   // 随机种子
};

// 工作线程函数
void worker_thread(queue<TestTask>& task_queue, atomic<int>& completed_tasks, 
                   int total_tasks, ofstream& out) {
    while (true) {
        // 从任务队列获取任务
        TestTask task;
        {
            lock_guard<mutex> lock(g_mutex);
            if (task_queue.empty()) {
                break;
            }
            task = task_queue.front();
            task_queue.pop();
        }
        
        int n = task.n;
        int round = task.round;
        uint32_t seed = task.seed;
        int mdeg = 3;
        int s = 1;

        // 生成图
        Graph g = rand_graph(n, mdeg, seed);

        // 计算参数
        auto lg = [](double x){ 
            return std::log(std::max(2.0, x)); 
        };

        double L = lg((double)n);
        int Sigma = std::max(1, (int)std::floor(std::pow(L, 1.0 / 3.0)));
        int Tau = std::max(1, (int)std::floor(std::pow(L, 2.0 / 3.0)));
        int ellTop = std::max(1, (int)std::ceil(L / std::max(1, Tau)));
        
        // 运行算法并计时
        auto t1 = std::chrono::steady_clock::now();
        BMSSPAlgo bm(g, s, g.n, Sigma, Tau);
        BMSSPResult res;
        bool bm_ok = true;
        std::string bm_err;
        try {
            res = bm.run(ellTop, INF_DIST, std::vector<int>{s});
        } catch (const std::exception &ex) {
            bm_ok = false;
            bm_err = ex.what();
        }
        auto t2 = std::chrono::steady_clock::now();

        double tb = std::chrono::duration<double>(t2 - t1).count();

        // 输出结果到文件
        {
            lock_guard<mutex> lock(g_mutex);
            out << fixed << setprecision(6);
            out << n/10000 << " " << tb*10000 << endl;
            out.flush();  // 确保及时写入
            
            // 更新进度
            completed_tasks++;
            int progress = (completed_tasks * 100) / total_tasks;
            cout << "Thread " << this_thread::get_id() 
                 << ": finished round " << round 
                 << " (n=" << n << ") in " << tb << "s. "
                 << "Progress: " << progress << "%" << endl;
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cout.tie(0);
    
    // 设置随机种子
    mt19937 rng(time(NULL));
    
    // 打开输出文件
    ofstream out("time.csv", ios::out | ios::app);
    out.tie(0);
    out << "x,y\n";
    
    // 创建任务队列
    queue<TestTask> task_queue;
    int total_tasks = 0;
    
    // 生成任务：从10000到10000000，步长10000
    for (int round = 10000; round <= 10000 * 1000; round += 10000) {
        int n = round;
        uint32_t seed = rng();  // 使用线程安全的随机数生成器
        
        task_queue.push({n, round, seed});
        total_tasks++;
    }
    
    cout << "Total tasks: " << total_tasks << endl;
    
    // 确定线程数（根据你的服务器核心数调整）
    // 获取硬件支持的并发线程数
    unsigned int num_threads = thread::hardware_concurrency();
    
    // 如果获取失败，使用默认值（建议使用物理核心数或略少于总逻辑核心数）
    if (num_threads == 0) {
        num_threads = 16;  // 设置为16个线程，略少于你的20核心
        cout << "Cannot determine hardware concurrency. Using " 
             << num_threads << " threads." << endl;
    } else {
        // 使用所有逻辑核心，或者减2留给系统使用
        num_threads = min(num_threads, (unsigned int)20);
        cout << "Detected " << num_threads << " concurrent threads available." << endl;
    }
    
    // 启动工作线程
    vector<thread> threads;
    atomic<int> completed_tasks(0);
    
    auto start_time = chrono::steady_clock::now();
    
    for (unsigned int i = 0; i < num_threads; i++) {
        threads.emplace_back(worker_thread, ref(task_queue), ref(completed_tasks), 
                            total_tasks, ref(out));
    }
    
    // 等待所有线程完成
    for (auto& t : threads) {
        t.join();
    }
    
    auto end_time = chrono::steady_clock::now();
    double total_time = chrono::duration<double>(end_time - start_time).count();
    
    cout << "\nAll tasks completed!" << endl;
    cout << "Total time: " << total_time << " seconds" << endl;
    cout << "Average time per task: " << total_time / total_tasks << " seconds" << endl;
    
    out.close();
    return 0;
}
