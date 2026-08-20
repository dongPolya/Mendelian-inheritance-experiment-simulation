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
            gene_info: 基因信息字典
            gene_locations: 基因位置字典
            lethal_genes: 致死配子字典，格式：
                {
                    "male": [("A",), ("B",)],  # 父本致死配子
                    "female": [("c",)]          # 母本致死配子
                }
            lethal_genotypes: 致死基因型列表（个体致死），如 ["AA", "Bb"]
                只要基因型包含这些组合就致死
        """
        self.name = name
        self.gene_info = gene_info if gene_info is not None else {}
        self.gene_locations = gene_locations if gene_locations is not None else {}

        # 初始化致死配子
        if lethal_genes is None:
            self.lethal_genes = {"male": [], "female": []}
        else:
            self.lethal_genes = lethal_genes

        # 初始化致死基因型
        self.lethal_genotypes = lethal_genotypes if lethal_genotypes is not None else []

    def get_trait(self, gene_order, genotype):
        """
        根据基因序号和基因型获取性状和显隐性
        """
        if gene_order not in self.gene_info:
            return None

        gene_data = self.gene_info[gene_order]
        genotype_trait_map = gene_data.get("genotype_trait_map", {})

        # 标准化基因型（排序，使 Aa 和 aA 视为相同）
        # 注意：不使用 reverse=True，以保证大写字母在前（与字典键一致）
        normalized_genotype = tuple(sorted(genotype))

        return genotype_trait_map.get(normalized_genotype)

    def get_phenotype_from_genotype(self, genotype_list):
        """
        根据完整基因型列表获取表型列表
        """
        phenotypes = []

        for idx, genotype in enumerate(genotype_list, start=1):
            trait_info = self.get_trait(idx, genotype)
            if trait_info and trait_info.get("trait"):
                phenotypes.append(trait_info["trait"])
            else:
                phenotypes.append("未知")

        return phenotypes

    def is_lethal_genotype(self, genotype_list):
        """
        判断基因型是否导致个体致死

        参数:
            genotype_list: 基因型列表，如 [("A", "A"), ("B", "b")]

        返回:
            bool: 是否致死
        """
        if not self.lethal_genotypes:
            return False

        # 将基因型展平为字符串
        flat_genotype = "".join("".join(gt) for gt in genotype_list)

        # 检查是否包含任何致死基因型
        for lethal_pattern in self.lethal_genotypes:
            if lethal_pattern in flat_genotype:
                return True

        return False

    def is_lethal_gamete(self, gamete, sex):
        """
        判断配子是否致死

        参数:
            gamete: 配子元组，如 ("A", "B", "c")
            sex: 性别，"male" 或 "female"

        返回:
            bool: 配子是否致死
        """
        lethal_list = self.lethal_genes.get(sex, [])

        for lethal_pattern in lethal_list:
            # 检查配子是否包含致死模式
            if self._gamete_matches_pattern(gamete, lethal_pattern):
                return True

        return False

    def _gamete_matches_pattern(self, gamete, pattern):
        """
        检查配子是否匹配致死模式

        参数:
            gamete: 配子元组，如 ("A", "B", "C")
            pattern: 致死模式元组，如 ("A", "B")

        返回:
            bool: 是否匹配
        """
        # 检查 pattern 中的所有元素是否都在 gamete 中
        return all(allele in gamete for allele in pattern)

    def get_location(self, gene_order, genotype):
        """
        获取基因在染色体上的位置
        """
        # 统一排序，不使用 reverse=True
        genotype_str = "".join(sorted(genotype))
        return self.gene_locations.get((gene_order, genotype_str), self.AUTOSOME)

    def cross_single_gene(self, parent1_genotype, parent2_genotype):
        """
        单基因杂交，计算后代基因型及比例
        """
        gametes1 = self._get_gametes(parent1_genotype)
        gametes2 = self._get_gametes(parent2_genotype)

        offspring = {}
        total = len(gametes1) * len(gametes2)

        for g1 in gametes1:
            for g2 in gametes2:
                # 标准化基因型，不使用 reverse=True
                genotype = tuple(sorted([g1, g2]))
                offspring[genotype] = offspring.get(genotype, 0) + 1

        for genotype in offspring:
            offspring[genotype] /= total

        return offspring

    def _get_gametes(self, genotype):
        """
        获取个体产生的配子
        """
        return list(genotype)

    def self_cross(self, individual):
        """
        自交：一个个体的自花授粉或自体受精
        """
        return self._cross(individual, individual, is_self=True)

    def hybrid_cross(self, individual1, individual2):
        """
        杂交：两个不同个体之间的交配
        """
        return self._cross(individual1, individual2, is_self=False)

    def free_mating(self, individuals):
        """
        自由交配：多个个体之间随机交配
        """
        all_offspring = []

        for i in range(len(individuals)):
            for j in range(i, len(individuals)):
                is_self = (i == j)
                offspring = self._cross(individuals[i], individuals[j], is_self=is_self)
                all_offspring.extend(offspring)

        return all_offspring

    def _cross(self, parent1, parent2, is_self=False):
        """
        内部方法：执行两个个体的杂交，考虑配子致死和个体致死

        参数:
            parent1: Individual 对象（父本）
            parent2: Individual 对象（母本）
            is_self: 是否为自交

        返回:
            Individual 对象列表：后代个体列表
        """
        if not parent1.genetic_type or not parent2.genetic_type:
            return []

        gene_count = len(parent1.genetic_type)

        # 生成所有可能的配子组合（考虑配子致死）
        parent1_gametes = self._generate_all_gametes(parent1.genetic_type)
        parent2_gametes = self._generate_all_gametes(parent2.genetic_type)

        # 过滤致死配子
        parent1_viable_gametes = [
            gamete for gamete in parent1_gametes
            if not self.is_lethal_gamete(gamete, "male")
        ]
        parent2_viable_gametes = [
            gamete for gamete in parent2_gametes
            if not self.is_lethal_gamete(gamete, "female")
        ]

        # 如果所有配子都致死，返回空列表
        if not parent1_viable_gametes or not parent2_viable_gametes:
            return []

        # 计算存活配子的概率（归一化）
        total_p1 = len(parent1_viable_gametes)
        total_p2 = len(parent2_viable_gametes)

        # 组合所有基因的后代
        combined_offspring = {}

        for gamete1 in parent1_viable_gametes:
            for gamete2 in parent2_viable_gametes:
                # 组合成基因型
                genotype_list = []
                for i in range(gene_count):
                    allele1 = gamete1[i]
                    allele2 = gamete2[i]
                    # 标准化基因型，不使用 reverse=True
                    genotype = tuple(sorted([allele1, allele2]))
                    genotype_list.append(genotype)

                genotype_tuple = tuple(genotype_list)

                # 计算该组合的概率
                probability = (1.0 / total_p1) * (1.0 / total_p2)

                combined_offspring[genotype_tuple] = (
                        combined_offspring.get(genotype_tuple, 0) + probability
                )

        # 创建后代个体，过滤致死个体
        offspring_list = []
        for genotype_list, probability in combined_offspring.items():
            # 检查个体是否致死
            if self.is_lethal_genotype(list(genotype_list)):
                continue

            # 获取表型
            phenotype = self.get_phenotype_from_genotype(list(genotype_list))

            # 创建后代个体
            offspring = Individual(
                genetic_type=[list(gt) for gt in genotype_list],
                phenotype=phenotype
            )
            offspring.probability = probability
            offspring_list.append(offspring)

        # 重新计算比例（去除致死个体后归一化）
        total_prob = sum(ind.probability for ind in offspring_list)
        if total_prob > 0:
            for ind in offspring_list:
                ind.probability /= total_prob

        return offspring_list

    def _generate_all_gametes(self, genotype_list):
        """
        生成所有可能的配子组合

        参数:
            genotype_list: 基因型列表，如 [("A", "a"), ("B", "b")]

        返回:
            list: 所有可能的配子列表，如 [("A", "B"), ("A", "b"), ("a", "B"), ("a", "b")]
        """
        if not genotype_list:
            return []

        # 从第一对基因开始
        result = [(allele,) for allele in genotype_list[0]]

        # 逐步添加后续基因的等位基因
        for i in range(1, len(genotype_list)):
            new_result = []
            for existing_gamete in result:
                for allele in genotype_list[i]:
                    new_gamete = existing_gamete + (allele,)
                    new_result.append(new_gamete)
            result = new_result

        return result

    def calculate_statistics(self, offspring_list):
        """
        计算后代的基因型和表型统计
        """
        genotype_ratio = {}
        phenotype_ratio = {}

        for individual in offspring_list:
            genotype_str = self._genotype_to_string(individual.genetic_type)
            genotype_ratio[genotype_str] = genotype_ratio.get(genotype_str, 0) + individual.probability

            phenotype_str = " | ".join(individual.phenotype) if individual.phenotype else "未知"
            phenotype_ratio[phenotype_str] = phenotype_ratio.get(phenotype_str, 0) + individual.probability

        return {
            "genotype_ratio": genotype_ratio,
            "phenotype_ratio": phenotype_ratio
        }

    def _genotype_to_string(self, genotype_list):
        """
        将基因型列表转换为字符串
        """
        return "".join("".join(gt) for gt in genotype_list)


class Individual:
    """个体类 - 表示具有特定基因型和表型的生物个体"""

    def __init__(self, genetic_type=None, phenotype=None, sex=None,
                 probability=None):
        """
        初始化个体
        """
        self.genetic_type = genetic_type if genetic_type is not None else []
        self.phenotype = phenotype if phenotype is not None else []
        self.sex = sex
        self.probability = probability if probability is not None else 1.0

    @property
    def genotype_ratio(self):
        """
        获取基因型占比
        """
        if not self.genetic_type:
            return {}

        genotype_str = "".join("".join(gt) for gt in self.genetic_type)
        return {genotype_str: self.probability}

    @property
    def phenotype_ratio(self):
        """
        获取表型占比
        """
        if not self.phenotype:
            return {}

        phenotype_str = " | ".join(self.phenotype)
        return {phenotype_str: self.probability}

    def get_info(self):
        """
        获取个体的详细信息
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
    # 创建遗传图谱：两对基因
    # 修复：补充了 ("A", "A") 的定义，删除了重复的 ("A", "a")
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
        }
    }

    gene_locations = {
        (1, "Aa"): GeneticMap.AUTOSOME,
        (2, "Bb"): GeneticMap.AUTOSOME
    }

    print("=" * 60)
    print("示例1：无致死情况的自交 - AaBb 自交")
    print("=" * 60)

    genetic_map_normal = GeneticMap(
        name="植物遗传图谱（无致死）",
        gene_info=gene_info,
        gene_locations=gene_locations
    )

    parent = Individual(
        genetic_type=[("A", "a"), ("B", "b")],
        phenotype=["宽叶", "粉叶"]
    )

    print(f"\n亲本信息:")
    print(parent.get_info())

    offspring_normal = genetic_map_normal.self_cross(parent)
    stats_normal = genetic_map_normal.calculate_statistics(offspring_normal)

    print(f"\n后代数量: {len(offspring_normal)}")
    print("\n表型比例:")
    for phenotype, ratio in sorted(stats_normal["phenotype_ratio"].items()):
        print(f"  {phenotype}: {ratio:.4f} ({ratio * 100:.2f}%)")

    print("\n" + "=" * 60)
    print("示例2：配子致死 - 父本A配子致死")
    print("=" * 60)

    lethal_genes_example = {
        "male": [("A",)],  # 父本的A配子致死
        "female": []  # 母本无致死配子
    }

    genetic_map_lethal_gamete = GeneticMap(
        name="植物遗传图谱（配子致死）",
        gene_info=gene_info,
        gene_locations=gene_locations,
        lethal_genes=lethal_genes_example
    )

    parent2 = Individual(
        genetic_type=[("A", "a"), ("B", "b")],
        phenotype=["宽叶", "粉叶"]
    )

    print(f"\n亲本信息:")
    print(parent2.get_info())
    print(f"\n致死配子设置: 父本A配子致死")

    offspring_lethal = genetic_map_lethal_gamete.self_cross(parent2)
    stats_lethal = genetic_map_lethal_gamete.calculate_statistics(offspring_lethal)

    print(f"\n后代数量: {len(offspring_lethal)}")
    print("\n表型比例:")
    for phenotype, ratio in sorted(stats_lethal["phenotype_ratio"].items()):
        print(f"  {phenotype}: {ratio:.4f} ({ratio * 100:.2f}%)")
    print("基因型比例:")
    for genotype, ratio in sorted(stats_lethal["genotype_ratio"].items()):
        print(f"  {genotype}: {ratio:.4f} ({ratio * 100:.2f}%)")

    print("\n" + "=" * 60)
    print("示例3：个体致死 - AA基因型致死")
    print("=" * 60)

    genetic_map_lethal_individual = GeneticMap(
        name="植物遗传图谱（个体致死）",
        gene_info=gene_info,
        gene_locations=gene_locations,
        lethal_genotypes=["AA"]  # 含有AA的个体致死
    )

    parent3 = Individual(
        genetic_type=[("A", "a"), ("B", "b")],
        phenotype=["宽叶", "粉叶"]
    )

    print(f"\n亲本信息:")
    print(parent3.get_info())
    print(f"\n致死基因型设置: AA致死")

    offspring_ind_lethal = genetic_map_lethal_individual.self_cross(parent3)
    stats_ind_lethal = genetic_map_lethal_individual.calculate_statistics(offspring_ind_lethal)

    print(f"\n后代数量: {len(offspring_ind_lethal)}")
    print("\n表型比例:")
    for phenotype, ratio in sorted(stats_ind_lethal["phenotype_ratio"].items()):
        print(f"  {phenotype}: {ratio:.4f} ({ratio * 100:.2f}%)")

    print("\n" + "=" * 60)
    print("示例4：同时存在配子致死和个体致死")
    print("=" * 60)

    lethal_genes_both = {
        "male": [("A",)],  # 父本A配子致死
        "female": [("b",)]  # 母本b配子致死
    }

    genetic_map_both = GeneticMap(
        name="植物遗传图谱（双重致死）",
        gene_info=gene_info,
        gene_locations=gene_locations,
        lethal_genes=lethal_genes_both,
        lethal_genotypes=["BB"]  # BB个体也致死
    )

    parent4 = Individual(
        genetic_type=[("A", "a"), ("B", "b")],
        phenotype=["宽叶", "粉叶"]
    )

    print(f"\n亲本信息:")
    print(parent4.get_info())
    print(f"\n致死设置:")
    print(f"  - 配子致死: 父本A配子，母本b配子")
    print(f"  - 个体致死: BB基因型")

    offspring_both = genetic_map_both.self_cross(parent4)

    if offspring_both:
        stats_both = genetic_map_both.calculate_statistics(offspring_both)

        print(f"\n后代数量: {len(offspring_both)}")
        print("\n表型比例:")
        for phenotype, ratio in sorted(stats_both["phenotype_ratio"].items()):
            print(f"  {phenotype}: {ratio:.4f} ({ratio * 100:.2f}%)")
    else:
        print("\n所有后代均致死！")

    print("\n" + "=" * 60)
    print("示例5：杂交实验 - 配子致死的影响")
    print("=" * 60)

    parent_male = Individual(
        genetic_type=[("A", "a"), ("B", "b")],
        phenotype=["宽叶", "粉叶"],
        sex="male"
    )

    parent_female = Individual(
        genetic_type=[("A", "a"), ("B", "b")],
        phenotype=["宽叶", "粉叶"],
        sex="female"
    )

    lethal_genes_hybrid = {
        "male": [("A",)],  # 父本A配子致死
        "female": []
    }

    genetic_map_hybrid = GeneticMap(
        name="杂交实验（配子致死）",
        gene_info=gene_info,
        gene_locations=gene_locations,
        lethal_genes=lethal_genes_hybrid
    )

    print(f"\n父本: {parent_male}")
    print(f"母本: {parent_female}")
    print(f"致死设置: 父本A配子致死")

    f1_offspring = genetic_map_hybrid.hybrid_cross(parent_male, parent_female)

    if f1_offspring:
        f1_stats = genetic_map_hybrid.calculate_statistics(f1_offspring)

        print(f"\nF1代数量: {len(f1_offspring)}")
        print("\nF1代表型比例:")
        for phenotype, ratio in sorted(f1_stats["phenotype_ratio"].items()):
            print(f"  {phenotype}: {ratio:.4f} ({ratio * 100:.2f}%)")
    else:
        print("\n所有F1代均致死！")
