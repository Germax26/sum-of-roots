import string
import sys
from typing import Dict, Generator, List, Tuple

DEGREE = int(sys.argv[1])

all_roots = string.ascii_letters
assert 0 <= DEGREE <= len(all_roots)

class Expr:
	def __init__(self, *sigmas: 'Expr'):
		self.sigmas: List[Sigma] = []
		for sigma in sigmas:
			if isinstance(sigma, Sigma): self.sigmas.append(sigma)
			elif isinstance(sigma, Expr): self.sigmas += sigma.sigmas
			else: assert False, "Unreachable"

	def __mul__(self, other: 'Expr | int'):
		if isinstance(other, Expr):
			new_sigmas = Expr()
			for a in self.sigmas:
				for b in other.sigmas:
					new_sigmas += a * b
			return Expr.collect(new_sigmas)
		elif isinstance(other, Sigma):
			return self * Expr(other)
		elif isinstance(other, int):
			return Expr(*[sigma * other for sigma in self.sigmas])

	def __add__(self, other: 'Expr') -> 'Expr': return Expr.collect(Expr(*self.sigmas, *other.sigmas))
	def __neg__(self): return self*-1
	def __sub__(self, other: 'Expr') -> 'Expr': return self+-other

	def __pow__(self, exp:int) -> 'Expr': 
		if exp == 0: return Sigma()
		assert exp > 0
		expr = self
		for i in range(exp-1): expr *= self
		return expr

	def __repr__(self):
		return " + ".join(map(str, self.sigmas))

	@staticmethod
	def collect(expr: 'Expr') -> 'Expr':
		forms: Dict[str, int] = {}
		for sigma in expr.sigmas: forms[sigma.form] = forms.get(sigma.form, 0) + sigma.coefficient
		return Expr(*[Sigma(form) * forms[form] for form in forms])

class Roots:
	def __init__(self, roots:str="", dict:Dict[str,int]={}):
		self.exps: Dict[str, int] = {}
		if dict:
			self.exps = dict
			return
		while roots:
			root = roots[0]
			roots = roots[1:]
			num = ""
			while roots and roots[0].isnumeric():
				num += roots[0]
				roots = roots[1:]
			self.exps[root] = self.exps.get(root, 0) + (int(num) if num else 1)

	def __getitem__(self, key: str) -> int:
		return self.exps[key]

	def get(self, key: str, default: int) -> int:
		return self.exps.get(key, default)

	def __mul__(self, other: 'Roots') -> 'Roots':
		return Roots(dict={root: self.get(root, 0) + other.get(root, 0) for root in self.exps | other.exps})

	def __pow__(self, other: int) -> 'Roots':
		return Roots(dict={root: self[root]*other for root in self.exps})

	def __repr__(self):
		return self.get_form()

	def get_reduced_form(self) -> str:
		exps: Dict[int, int] = {}
		for root in self.exps: exps[self.exps[root]] = exps.get(self.exps[root], 0) + 1

		form = ""
		i = 0
		for exp in reversed(sorted(exps)):
			for j in range(exps[exp]):
				form += all_roots[i] + (str(exp) if exp != 1 else "")
				i += 1

		return form

	def get_form(self) -> str:
		return "".join([root + (str(self.exps[root]) if self.exps[root] != 1 else "") for root in self.exps])
		
	@staticmethod
	def get_from(form: str="") -> Generator['Roots', None, None]:
		if not form: 
			yield Roots()
			return

		roots_exps = Roots(form).exps

		exps: Dict[int, int] = {}
		for root in roots_exps: exps[roots_exps[root]] = exps.get(roots_exps[root], 0) + 1
		
		def choose(objects: str, k:int) -> Generator[str, None, None]:
			if k:
				for i in range(len(objects)):
					for combo in choose(objects[i+1:],k-1): yield objects[i]+combo
			else: yield ""

		def foo(objects: str, ks: List[Tuple[int, int]]) -> Generator[Roots, None, None]:
			if not len(ks): 
				yield Roots()
				return
			exp, k, ksr = ks[0][0], ks[0][1], ks[1:]
			sobjects = set(objects)
			for choice in choose(objects, k):
				new_objs = "".join(list(sobjects - set(choice)))
				for subchoice in foo(new_objs, ksr):
					yield Roots("".join(choice)) ** exp * subchoice

		yield from foo(all_roots[:DEGREE], list(exps.items()))

class Sigma(Expr):
	def __init__(self, form="", coefficient=1):
		self.coefficient: int = coefficient
		self.form: str = form
		super().__init__(self)

	def __repr__(self):
		if self.coefficient == 0: return "0"
		coefficient = {1:"",-1:"-"}.get(self.coefficient, str(self.coefficient))
		return coefficient + "Σ" + self.form

	def __mul__(self, other: Expr | int) -> Expr:
		if isinstance(other, int): return Sigma(self.form, self.coefficient*other)
		if not isinstance(other, Sigma): return super().__mul__(other)
		
		coeff = self.coefficient * other.coefficient
		if coeff == 0: return Sigma(coefficient=0)
		roots: Dict[str, int] = {}

		# TODO: This function needs to be optimised.
		# Currently, it is generating all the possible root combinations, multiplyin them and then collecting like terms.
		# I've looked at the expansions of the (Σa)^N series, and there seems to be a connection between the coefficients and combinatorics
		# There should be a combinatoric way to multiple two of these Sigma's together.
		# But as of now, I am not sure what it is.

		for root_self in Roots.get_from(self.form):
			for root_other in Roots.get_from(other.form):
				root = root_self * root_other
				form = root.get_reduced_form()
				roots[form] = roots.get(form, 0) + 1

		for form in roots:
			all_from = [x for x in Roots.get_from(form)]
			amount = roots[form] / len(all_from)
			int_amount = int(amount)
			assert int_amount == amount, (int_amount, amount) # symmetry of roots. I have not proved that this will always be true, but I hope it is
			roots[form] = int_amount

		return Expr(*[Sigma(form, coeff*roots[form]) for form in roots])
