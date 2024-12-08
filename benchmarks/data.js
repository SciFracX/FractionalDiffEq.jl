window.BENCHMARK_DATA = {
  "lastUpdate": 1733672375965,
  "repoUrl": "https://github.com/SciFracX/FractionalDiffEq.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "48be1127311bdb4798db4711a0e422cc68644703",
          "message": "Merge pull request #108 from ErikQQY/qqy/better_bench\n\nAdd better benchmarks",
          "timestamp": "2024-04-06T00:50:44+08:00",
          "tree_id": "92b826bff9fa651e84c7a88f3ce70618e419c553",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/48be1127311bdb4798db4711a0e422cc68644703"
        },
        "date": 1712335993264,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 18587804,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 18020466,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419536\nallocs=264997\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 18078300,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22337024\nallocs=264985\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "37963e0f168edc14574c98be212fcc83bdd31355",
          "message": "Merge pull request #109 from ErikQQY/qqy/triger_bench\n\nTrigger benchmarks",
          "timestamp": "2024-04-06T12:58:39+08:00",
          "tree_id": "808e5a17757f6632791169aa8d1c15bc9683f7c7",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/37963e0f168edc14574c98be212fcc83bdd31355"
        },
        "date": 1712379651197,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17310408.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16730695.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419536\nallocs=264997\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16781819.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22337024\nallocs=264985\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "302ce2a3de9d5c0f72d40c2ee1d4299bc55c1a17",
          "message": "Merge pull request #110 from ErikQQY/qqy/better_ben\n\nTest benchmarks",
          "timestamp": "2024-04-06T13:11:05+08:00",
          "tree_id": "5e438e322d22539c126c00d7af0c32a3bb3a264b",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/302ce2a3de9d5c0f72d40c2ee1d4299bc55c1a17"
        },
        "date": 1712380402650,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 18292251,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 17864163.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419536\nallocs=264997\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 17813185.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a101dbf65fb6651424cec395f5861ed45d61fbed",
          "message": "Merge pull request #111 from ErikQQY/qqy/test_bench\n\nTest CI benchmarks",
          "timestamp": "2024-04-06T13:30:43+08:00",
          "tree_id": "0967cfb08441a2e61a4a0cb00faa22b8423f7403",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/a101dbf65fb6651424cec395f5861ed45d61fbed"
        },
        "date": 1712381578766,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17007258.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16916692,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419184\nallocs=264978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16571171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "61a089ca536331e7dfa346e1d3554aa8474da65d",
          "message": "Merge pull request #112 from ErikQQY/qqy/test_bench\n\nTrigger again",
          "timestamp": "2024-06-01T20:39:28+08:00",
          "tree_id": "ebf3007828d12bd460e69e68a28e29deab339737",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/61a089ca536331e7dfa346e1d3554aa8474da65d"
        },
        "date": 1717245718580,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17357485,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22961808\nallocs=265810\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16803226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22417184\nallocs=264853\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16816974.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22334672\nallocs=264841\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "4c10b03f697e027bf99839b130a372a730392d2e",
          "message": "Add formatter\n\nSigned-off-by: ErikQQY <2283984853@qq.com>",
          "timestamp": "2024-06-01T20:41:21+08:00",
          "tree_id": "89bf8adc6516b1ca4ee251500427894761c01381",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/4c10b03f697e027bf99839b130a372a730392d2e"
        },
        "date": 1717245886944,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17545076,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16847063.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419184\nallocs=264978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16937641.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "d245790c5f3b26d675ddca05585d84364f4454f4",
          "message": "format\n\nSigned-off-by: ErikQQY <2283984853@qq.com>",
          "timestamp": "2024-06-01T21:31:31+08:00",
          "tree_id": "795b6e4e084210be5326d56d92ae6c758cbdc5f4",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/d245790c5f3b26d675ddca05585d84364f4454f4"
        },
        "date": 1717248832892,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17002887.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16275364,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419184\nallocs=264978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16414135,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "6a756727a97441ffff86487db1165b6c7e02d787",
          "message": "Update README\n\nSigned-off-by: ErikQQY <2283984853@qq.com>",
          "timestamp": "2024-06-01T21:36:42+08:00",
          "tree_id": "fae136515a245054fa88d11492de8578d5c9bd75",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/6a756727a97441ffff86487db1165b6c7e02d787"
        },
        "date": 1717249145010,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 18096207,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 17524745.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419184\nallocs=264978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 17460771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "973c3db818b739e7f2a2327a5f3e16efd846626b",
          "message": "Update docs\n\nSigned-off-by: ErikQQY <2283984853@qq.com>",
          "timestamp": "2024-06-25T23:04:09+08:00",
          "tree_id": "6caf46537a1d4b0df4c1a3ccceeadee5e7b7538e",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/973c3db818b739e7f2a2327a5f3e16efd846626b"
        },
        "date": 1719328005479,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 20946614,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22961808\nallocs=265810\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16455634,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22417184\nallocs=264853\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16555383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22334672\nallocs=264841\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "8576eb05cb1aa9b3d30d3becaa2b473f3435b091",
          "message": "Update FDDEProblem\n\nSigned-off-by: ErikQQY <2283984853@qq.com>",
          "timestamp": "2024-06-26T11:09:21+08:00",
          "tree_id": "7e8322a04d9d09f2ec410df574b39dd25a10edc5",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/8576eb05cb1aa9b3d30d3becaa2b473f3435b091"
        },
        "date": 1719371505375,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 20342051,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16888141,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419184\nallocs=264978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16795517,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "eccc0254d0199669b5a0282b5a79162098a5b498",
          "message": "Merge pull request #114 from alexfikl/fix-brusselator\n\nFix Brusselator right-hand side in some examples",
          "timestamp": "2024-08-23T00:37:44+08:00",
          "tree_id": "7d9841baa5926256c904fddf4dcb4fb72343cf8c",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/eccc0254d0199669b5a0282b5a79162098a5b498"
        },
        "date": 1724344850305,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 19841908,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22206432\nallocs=252355\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 18515591.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21712560\nallocs=252308\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 18477555,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21649600\nallocs=252648\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "904a32232199f8ce9a224525f084f29fe37bc950",
          "message": "Merge pull request #115 from ErikQQY/qqy/drop_unpack\n\nDrop UnPack",
          "timestamp": "2024-09-19T00:03:03+08:00",
          "tree_id": "104fb32c9e9682a061737f7106493e3582ca8b06",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/904a32232199f8ce9a224525f084f29fe37bc950"
        },
        "date": 1726675526663,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 19422373,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22204432\nallocs=252230\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 18519249,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21710560\nallocs=252183\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 18570106,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21647600\nallocs=252523\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "0bfce8b265049f1fea81336afb3b3ced92f495c9",
          "message": "update README\n\nSigned-off-by: Qingyu Qu <2283984853@qq.com>",
          "timestamp": "2024-09-19T16:23:23+08:00",
          "tree_id": "d95d641633915b6c686c89df91c58650da62589b",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/0bfce8b265049f1fea81336afb3b3ced92f495c9"
        },
        "date": 1726734666938,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 19542759.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22204432\nallocs=252230\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 18801598,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21710560\nallocs=252183\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 18617032,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21647600\nallocs=252523\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9a7b5d4b0972bb896779c20c65b1ad6745808ac4",
          "message": "Merge pull request #117 from ErikQQY/qqy/up\n\nUpdate compat",
          "timestamp": "2024-10-09T01:31:35+08:00",
          "tree_id": "d9bf76be02b24a72c4727745b57560d8379dc97f",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/9a7b5d4b0972bb896779c20c65b1ad6745808ac4"
        },
        "date": 1728408873431,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 18397976.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23069952\nallocs=385981\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 17944961,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22575856\nallocs=385922\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 17901242,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22512816\nallocs=386357\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "8dd594dc33ad7a54dc487a0ae24c32248f15eb3e",
          "message": "Merge pull request #118 from ErikQQY/qqy/initial_fode\n\nInitial improvement for FODE solvers",
          "timestamp": "2024-10-09T23:26:58+08:00",
          "tree_id": "50cdaaba2ecbd263a44e900b4f6059dd88a1d4bd",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/8dd594dc33ad7a54dc487a0ae24c32248f15eb3e"
        },
        "date": 1728487795381,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 16727289,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20414480\nallocs=338020\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16744188,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20559280\nallocs=340213\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16152978.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20262992\nallocs=318999\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "eb87249d910a600044b1d9020256e73a831b7729",
          "message": "Merge pull request #120 from ErikQQY/qqy/remove\n\nRemove redundant problems and solvers",
          "timestamp": "2024-10-26T17:05:13+08:00",
          "tree_id": "f7e78b23915bb13130a98a333f0bd63d505bbd08",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/eb87249d910a600044b1d9020256e73a831b7729"
        },
        "date": 1729933674175,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 16573801.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20414480\nallocs=338020\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16448950.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20559280\nallocs=340213\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 15710771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20262992\nallocs=318999\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9801d8e44365448c5eea8005ac98d7fe55d54f99",
          "message": "Merge pull request #121 from ErikQQY/qqy/ext_oop\n\nFix oop suport for extension",
          "timestamp": "2024-10-26T20:24:36+08:00",
          "tree_id": "eb527a2df87e2da3d32555958e15311673ebe998",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/9801d8e44365448c5eea8005ac98d7fe55d54f99"
        },
        "date": 1729945652914,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17354218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20414480\nallocs=338020\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 17163478.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20559280\nallocs=340213\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16402366,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20262992\nallocs=318999\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "e15d6e48ca9850c4e6ea39f3c8660e35fdb911ad",
          "message": "Fix multiterms PECE methods",
          "timestamp": "2024-12-08T23:36:42+08:00",
          "tree_id": "e855e168c62d18b1d433cad7f2e82078541fbe2f",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/e15d6e48ca9850c4e6ea39f3c8660e35fdb911ad"
        },
        "date": 1733672374397,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 15721938,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20309760\nallocs=335848\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 15644044,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20454560\nallocs=338041\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 14336648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19959280\nallocs=312751\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}