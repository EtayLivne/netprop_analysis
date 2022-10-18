# print("-1")
from pipelines.infection_resiliency_pipeline.run_pipeline import run_once
# from ray import remote, get
# from threading import Thread
# print(0)
#
# @remote
# def task_say(x):
#     print(x)
#
# def thread_say(x):
#     print("thread")
#     print(x)
#
# @remote
# class Actor:
#     def __init__(self, x):
#         self.x = x
#
#     def say(self):
#         #print(self.x)
#         #get(task_say.remote(self.x))
#         Thread(target=thread_say, args=(self.x,)).start()
#
# @remote
# class ActorFactory:
#     def get_actor(self, x):
#         return Actor.options(max_concurrency=2).remote(x)
#
#
# def exp():
#     factory = ActorFactory.remote()
#     actors = [get(factory.get_actor.remote(i)) for i in range(5)]
#     for a in actors:
#         get(a.say.remote())
#
# print("shoes")
#
#
#
#
def main() -> None:
    individual_interactors_metadata = "/data/metadata.json"
    prop_results_file = "/data/all.csv"
    resiliency_dir = "/output"
    # individual_interactors_metadata = "/data/propagations/krogan_interactors/individual_interactors/metadata.json"
    # prop_results_file = "/data/propagations/krogan_interactors/individual_interactors/all.csv"
    # resiliency_dir = "/data/resiliency"
    ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    min_subset_size = 5
    max_subset_size = 30
    top_propagated_rank_threshold = 1000
    top_resilient_rank_threshold = 100
    run_name = "test"
    print("TEST!")
    run_once(individual_interactors_metadata, prop_results_file, resiliency_dir,
             ratios, min_subset_size, max_subset_size, top_propagated_rank_threshold, top_resilient_rank_threshold,
             run_name, reset_state=True)
#
#
#
# print(1)
if __name__ == "__main__":
    print(2)
    main()

# print("yo")