from collections import Counter
import math


class Cluster():
    def __init__(self, item: tuple[int, str], threshold) -> None:
        self.map_index_member = {}
        self.map_index_member[item[0]] = item[1]
        self.threshold = threshold

    def __lt__(self, other):
        """
        Compare Cluster objects based on size of each object
        """
        return self.get_size() < other.get_size()

    def match_score(self, s1: str, s2: str) -> float:
        return float(sum(c1 == c2 for c1, c2 in zip(s1, s2)))/float(len(s1))

    def add(self, item: tuple[int, str]):
        """
        the item is added if the distances from it to all the members are lower than the threshold
        """
        idx, val = item
        for _, v in self.map_index_member.items():
            if self.match_score(v, val) <= self.threshold:
                return
        self.map_index_member[idx] = val
        return

    def get_indices(self):
        return list(self.map_index_member.keys())

    def get_size(self):
        return len(self.map_index_member)

    def get_members(self):
        return set(self.map_index_member.values())

    def get_representative(self):
        """
        Representative is the most common member in the cluster
        """
        return Counter(self.map_index_member.values()).most_common(1)[0][0]

    def print_cluster(self):
        """
        Print the members of cluster
        """
        out = ""
        for k, v in self.map_index_member.items():
            score = self.match_score(v, self.get_representative())
            if math.isclose(score, 1):
                score = "*"
            else:
                score = "{:.2f}".format(score)
            out += f"{k} {v} {score}\n"
        return out
