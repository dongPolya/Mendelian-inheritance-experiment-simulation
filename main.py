class GeneticMap:
    """遗传图谱类 - 管理遗传规则、基因关系和显隐性逻辑"""

    # 显隐性类型常量
    DOMINANT = "显性"
    RECESSIVE = "隐性"
    INCOMPLETE_DOMINANT = "不完全显性"

    # 染色体类型常量
    AUTOSOME = "常染色体"
    X_CHROMOSOME = "X染色体"
    Y_CHROMOSOME = "Y染色体"

    def __init__(self, name=None, gene_info=None, gene_locations=None,
                 lethal_genes=None, lethal_genotypes=None):
        """
        初始化遗传图谱

        参数:
            name: 遗传图谱名称
            gene_info: 基因信息字典，格式：
                {
                    1: {
                        "alleles": ["A", "a"],
                        "genotype_trait_map": {
                            ("A", "A"): {"dominance": "显性", "trait": "宽叶"},
                            ("A", "a"): {"dominance": "显性", "trait": "宽叶"},
                            ("a", "a"): {"dominance": "隐性", "trait": "窄叶"}
                        }
                    },
                    2: {
                        "alleles": ["B", "b"],
                        "genotype_trait_map": {
                            ("B", "B"): {"dominance": "显性", "trait": "红叶"},
                            ("B", "b"): {"dominance": "不完全显性", "trait": "粉叶"},
                            ("b", "b"): {"dominance": "隐性", "trait": "黄叶"}
                        }
                    }
                }
            gene_locations: 基因位置字典，格式：
                {
                    (1, "Aa"): "常染色体",
                    (2, "Bb"): "X染色体"
                }
            lethal_genes: 致死基因列表，如 ["AA", "bb"]
            lethal_genotypes: 致死基因型列表，如 [("A", "A", "B", "B")]
        """
        self.name = name
        self.gene_info = gene_info if gene_info is not None else {}
        self.gene_locations = gene_locations if gene_locations is not None else {}
        self.lethal_genes = lethal_genes if lethal_genes is not None else []
        self.lethal_genotypes = lethal_genotypes if lethal_genotypes is not None else []

    def get_trait(self, gene_order, genotype):
        """
        根据基因序号和基因型获取性状和显隐性

        参数:
            gene_order: 基因序号 (int)
            genotype: 基因型元组，如 ("A", "a")

        返回:
            dict: {"dominance": 显隐性类型, "trait": 性状名称} 或 None
        """
        if gene_order not in self.gene_info:
            return None

        gene_data = self.gene_info[gene_order]
        genotype_trait_map = gene_data.get("genotype_trait_map", {})

        # 标准化基因型（排序，使 Aa 和 aA 视为相同）
        normalized_genotype = tuple(sorted(genotype, reverse=True))

        return genotype_trait_map.get(normalized_genotype)

    def get_phenotype_from_genotype(self, genotype_list):
        """
        根据完整基因型列表获取表型列表

        参数:
            genotype_list: 基因型列表，如 [("A", "a"), ("B", "b")]

        返回:
            list: 表型列表，如 ["宽叶", "粉叶"]
        """
        phenotypes = []

        for idx, genotype in enumerate(genotype_list, start=1):
            trait_info = self.get_trait(idx, genotype)
            if trait_info and trait_info.get("trait"):
                phenotypes.append(trait_info["trait"])
            else:
                phenotypes.append("未知")

        return phenotypes

    def is_lethal(self, genotype_list):
        """
        判断基因型是否致死

        参数:
            genotype_list: 基因型列表

        返回:
            bool: 是否致死
        """
        # 检查单个基因是否致死
        for genotype in genotype_list:
            genotype_str = "".join(genotype)
            if genotype_str in self.lethal_genes:
                return True

        # 检查组合基因型是否致死
        flat_genotype = tuple()
        for genotype in genotype_list:
            flat_genotype += genotype

        if flat_genotype in self.lethal_genotypes:
            return True

        return False

    def get_location(self, gene_order, genotype):
        """
        获取基因在染色体上的位置

        参数:
            gene_order: 基因序号
            genotype: 基因型字符串，如 "Aa"

        返回:
            str: 染色体位置
        """
        genotype_str = "".join(sorted(genotype, reverse=True))
        return self.gene_locations.get((gene_order, genotype_str), self.AUTOSOME)

    def cross_single_gene(self, parent1_genotype, parent2_genotype):
        """
        单基因杂交，计算后代基因型及比例

        参数:
            parent1_genotype: 父本基因型元组，如 ("A", "a")
            parent2_genotype: 母本基因型元组，如 ("A", "a")

        返回:
            dict: {基因型元组: 比例}，如 {("A", "A"): 0.25, ("A", "a"): 0.5, ("a", "a"): 0.25}
        """
        # 获取配子
        gametes1 = self._get_gametes(parent1_genotype)
        gametes2 = self._get_gametes(parent2_genotype)

        # 计算后代基因型
        offspring = {}
        total = len(gametes1) * len(gametes2)

        for g1 in gametes1:
            for g2 in gametes2:
                # 标准化基因型
                genotype = tuple(sorted([g1, g2], reverse=True))
                offspring[genotype] = offspring.get(genotype, 0) + 1

        # 转换为比例
        for genotype in offspring:
            offspring[genotype] /= total

        return offspring

    def _get_gametes(self, genotype):
        """
        获取个体产生的配子

        参数:
            genotype: 基因型元组

        返回:
            list: 配子列表
        """
        return list(genotype)

    def self_cross(self, individual):
        """
        自交：一个个体的自花授粉或自体受精

        参数:
            individual: Individual 对象

        返回:
            Individual 对象列表：后代个体列表（已去除致死个体）
        """
        return self._cross(individual, individual)

    def hybrid_cross(self, individual1, individual2):
        """
        杂交：两个不同个体之间的交配

        参数:
            individual1: Individual 对象（父本）
            individual2: Individual 对象（母本）

        返回:
            Individual 对象列表：后代个体列表（已去除致死个体）
        """
        return self._cross(individual1, individual2)

    def free_mating(self, individuals):
        """
        自由交配：多个个体之间随机交配

        参数:
            individuals: Individual 对象列表

        返回:
            Individual 对象列表：所有可能的后代个体列表（已去除致死个体）
        """
        all_offspring = []

        # 所有个体两两交配（包括自交）
        for i in range(len(individuals)):
            for j in range(i, len(individuals)):
                offspring = self._cross(individuals[i], individuals[j])
                all_offspring.extend(offspring)

        return all_offspring

    def _cross(self, parent1, parent2):
        """
        内部方法：执行两个个体的杂交

        参数:
            parent1: Individual 对象
            parent2: Individual 对象

        返回:
            Individual 对象列表：后代个体列表
        """
        if not parent1.genetic_type or not parent2.genetic_type:
            return []

        # 对每对基因进行杂交
        gene_count = len(parent1.genetic_type)
        all_gene_offspring = []

        for i in range(gene_count):
            p1_genotype = parent1.genetic_type[i]
            p2_genotype = parent2.genetic_type[i]

            offspring_ratio = self.cross_single_gene(p1_genotype, p2_genotype)
            all_gene_offspring.append(offspring_ratio)

        # 组合所有基因的后代
        combined_offspring = self._combine_genes(all_gene_offspring)

        # 创建后代个体
        offspring_list = []
        for genotype_list, probability in combined_offspring.items():
            # 检查是否致死
            if self.is_lethal(genotype_list):
                continue

            # 获取表型
            phenotype = self.get_phenotype_from_genotype(genotype_list)

            # 创建后代个体
            offspring = Individual(
                genetic_type=list(genotype_list),
                phenotype=phenotype
            )
            offspring.probability = probability
            offspring_list.append(offspring)

        # 重新计算比例（去除致死个体后）
        total_prob = sum(ind.probability for ind in offspring_list)
        if total_prob > 0:
            for ind in offspring_list:
                ind.probability /= total_prob

        return offspring_list

    def _combine_genes(self, all_gene_offspring):
        """
        组合多对基因的后代

        参数:
            all_gene_offspring: 每对基因的后代比例列表

        返回:
            dict: {(基因型1, 基因型2, ...): 比例}
        """
        if not all_gene_offspring:
            return {}

        # 从第一对基因开始
        result = {tuple(): 1.0}

        for gene_offspring in all_gene_offspring:
            new_result = {}

            for existing_genotype, existing_prob in result.items():
                for new_genotype, new_prob in gene_offspring.items():
                    combined = existing_genotype + (new_genotype,)
                    new_result[combined] = new_result.get(combined, 0) + existing_prob * new_prob

            result = new_result

        return result

    def calculate_statistics(self, offspring_list):
        """
        计算后代的基因型和表型统计

        参数:
            offspring_list: Individual 对象列表

        返回:
            dict: {
                "genotype_ratio": {基因型字符串: 比例},
                "phenotype_ratio": {表型字符串: 比例}
            }
        """
        genotype_ratio = {}
        phenotype_ratio = {}

        for individual in offspring_list:
            # 基因型统计
            genotype_str = self._genotype_to_string(individual.genetic_type)
            genotype_ratio[genotype_str] = genotype_ratio.get(genotype_str, 0) + individual.probability

            # 表型统计
            phenotype_str = " | ".join(individual.phenotype) if individual.phenotype else "未知"
            phenotype_ratio[phenotype_str] = phenotype_ratio.get(phenotype_str, 0) + individual.probability

        return {
            "genotype_ratio": genotype_ratio,
            "phenotype_ratio": phenotype_ratio
        }

    def _genotype_to_string(self, genotype_list):
        """
        将基因型列表转换为字符串

        参数:
            genotype_list: 基因型列表，如 [("A", "a"), ("B", "b")]

        返回:
            str: 基因型字符串，如 "AaBb"
        """
        return "".join("".join(genotype) for genotype in genotype_list)


class Individual:
    """个体类 - 表示具有特定基因型和表型的生物个体"""

    def __init__(self, genetic_type=None, phenotype=None, sex=None,
                 probability=None):
        """
        初始化个体

        参数:
            genetic_type: 基因型列表，如 [("A", "a"), ("B", "b")]
            phenotype: 表型列表，如 ["宽叶", "粉叶"]
            sex: 性别，"male" 或 "female"
            probability: 该个体出现的概率（用于统计）
        """
        self.genetic_type = genetic_type if genetic_type is not None else []
        self.phenotype = phenotype if phenotype is not None else []
        self.sex = sex
        self.probability = probability if probability is not None else 1.0

    @property
    def genotype_ratio(self):
        """
        获取基因型占比（对于单个个体，就是其概率）

        返回:
            dict: {基因型字符串: 比例}
        """
        if not self.genetic_type:
            return {}

        genotype_str = "".join("".join(gt) for gt in self.genetic_type)
        return {genotype_str: self.probability}

    @property
    def phenotype_ratio(self):
        """
        获取表型占比（对于单个个体，就是其概率）

        返回:
            dict: {表型字符串: 比例}
        """
        if not self.phenotype:
            return {}

        phenotype_str = " | ".join(self.phenotype)
        return {phenotype_str: self.probability}

    def get_info(self):
        """
        获取个体的详细信息

        返回:
            str: 个体信息字符串
        """
        genotype_str = "".join("".join(gt) for gt in self.genetic_type) if self.genetic_type else "未知"
        phenotype_str = " | ".join(self.phenotype) if self.phenotype else "未知"
        sex_str = "雄性" if self.sex == "male" else "雌性" if self.sex == "female" else "未知"

        return (f"基因型: {genotype_str}\n"
                f"表型: {phenotype_str}\n"
                f"性别: {sex_str}\n"
                f"概率: {self.probability:.4f}")

    def __str__(self):
        """字符串表示"""
        genotype_str = "".join("".join(gt) for gt in self.genetic_type) if self.genetic_type else "?"
        phenotype_str = "/".join(self.phenotype) if self.phenotype else "?"
        return f"Individual({genotype_str}, {phenotype_str}, P={self.probability:.4f})"


# ==================== 使用示例 ====================

if __name__ == "__main__":
    # 创建遗传图谱：两对基因，一对完全显性，一对不完全显性
    gene_info = {
        1: {
            "alleles": ["A", "a"],
            "genotype_trait_map": {
                ("A", "A"): {"dominance": GeneticMap.DOMINANT, "trait": "宽叶"},
                ("A", "a"): {"dominance": GeneticMap.DOMINANT, "trait": "宽叶"},
                ("a", "a"): {"dominance": GeneticMap.RECESSIVE, "trait": "窄叶"}
            }
        },
        2: {
            "alleles": ["B", "b"],
            "genotype_trait_map": {
                ("B", "B"): {"dominance": GeneticMap.DOMINANT, "trait": "红叶"},
                ("B", "b"): {"dominance": GeneticMap.INCOMPLETE_DOMINANT, "trait": "粉叶"},
                ("b", "b"): {"dominance": GeneticMap.RECESSIVE, "trait": "黄叶"}
            }
        },
        3: {
            "alleles": ["C", "c"],
            "genotype_trait_map": {
                ("C", "C"): {"dominance": GeneticMap.DOMINANT, "trait": "高叶"},
                ("C", "c"): {"dominance": GeneticMap.DOMINANT, "trait": "高叶"},
                ("c", "c"): {"dominance": GeneticMap.RECESSIVE, "trait": "矮叶"}
            }
        }
    }

    gene_locations = {
        (1, "Aa"): GeneticMap.AUTOSOME,
        (2, "Bb"): GeneticMap.AUTOSOME,
        (3, "Cc"): GeneticMap.AUTOSOME
    }

    genetic_map = GeneticMap(
        name="植物遗传图谱",
        gene_info=gene_info,
        gene_locations=gene_locations,
    )

    print("=" * 60)
    print("示例1：自交实验 - AaBb 自交")
    print("=" * 60)

    # 创建一个基因型为 AaBb 的个体
    parent = Individual(
        genetic_type=[("A", "a"), ("B", "b"), ("C", "c")],
        phenotype=["宽叶", "粉叶", "高叶"],
        sex=None
    )

    print(f"\n亲本信息:")
    print(parent.get_info())

    # 自交
    offspring = genetic_map.self_cross(parent)

    print(f"\n后代数量: {len(offspring)}")
    print("\n后代详情:")
    for i, ind in enumerate(offspring, 1):
        print(f"\n后代 {i}:")
        print(ind.get_info())

    # 统计
    stats = genetic_map.calculate_statistics(offspring)

    print("\n" + "=" * 60)
    print("基因型比例:")
    for genotype, ratio in sorted(stats["genotype_ratio"].items()):
        print(f"  {genotype}: {ratio:.4f} ({ratio * 100:.2f}%)")

    print("\n表型比例:")
    for phenotype, ratio in sorted(stats["phenotype_ratio"].items()):
        print(f"  {phenotype}: {ratio:.4f} ({ratio * 100:.2f}%)")

    print("\n" + "=" * 60)
    print("示例2：杂交实验 - AABB × aabb")
    print("=" * 60)

    parent1 = Individual(
        genetic_type=[("A", "A"), ("B", "B")],
        phenotype=["宽叶", "红叶"]
    )

    parent2 = Individual(
        genetic_type=[("a", "a"), ("b", "b")],
        phenotype=["窄叶", "黄叶"]
    )

    print(f"\n父本: {parent1}")
    print(f"母本: {parent2}")

    f1_offspring = genetic_map.hybrid_cross(parent1, parent2)

    print(f"\nF1代:")
    for ind in f1_offspring:
        print(ind.get_info())

    print("\n" + "=" * 60)
    print("示例3：F1自交产生F2代")
    print("=" * 60)

    if f1_offspring:
        f2_offspring = genetic_map.self_cross(f1_offspring[0])

        print(f"F2代数量: {len(f2_offspring)}")

        f2_stats = genetic_map.calculate_statistics(f2_offspring)

        print("\nF2代表型比例:")
        for phenotype, ratio in sorted(f2_stats["phenotype_ratio"].items()):
            print(f"  {phenotype}: {ratio:.4f} ({ratio * 100:.2f}%)")

    print("\n" + "=" * 60)
    print("示例4：自由交配")
    print("=" * 60)

    # 创建多个个体
    individuals = [
        Individual(genetic_type=[("A", "A"), ("B", "B")], phenotype=["宽叶", "红叶"]),
        Individual(genetic_type=[("A", "a"), ("B", "b")], phenotype=["宽叶", "粉叶"]),
        Individual(genetic_type=[("a", "a"), ("b", "b")], phenotype=["窄叶", "黄叶"])
    ]

    free_offspring = genetic_map.free_mating(individuals)

    print(f"自由交配后代数量: {len(free_offspring)}")

    free_stats = genetic_map.calculate_statistics(free_offspring)

    print("\n自由交配后代表型比例:")
    for phenotype, ratio in sorted(free_stats["phenotype_ratio"].items()):
        print(f"  {phenotype}: {ratio:.4f} ({ratio * 100:.2f}%)")
