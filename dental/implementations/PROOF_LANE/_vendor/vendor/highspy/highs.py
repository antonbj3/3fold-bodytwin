from __future__ import annotations
import sys
import re
import numpy as np
from numbers import Integral
from itertools import product
from threading import Thread, local, RLock, Lock
from typing import Literal, Optional, Any, overload, Callable, Sequence, Mapping, Iterable, Iterator, SupportsIndex, cast, Union, Tuple, TYPE_CHECKING
from weakref import proxy
if sys.version_info >= (3, 10):
    from typing import TypeAlias
else:
    from typing_extensions import TypeAlias
from ._core import ObjSense, HighsVarType, HighsStatus, cb, _Highs, kHighsInf
if TYPE_CHECKING:
    ndarray_object_type = np.ndarray[Any, np.dtype[np.object_]]
    HighspyLeafTypes: TypeAlias = Union[int, Integral, 'highs_var', 'highs_cons', 'highs_linear_expression']
    HighspyNestedIndex: TypeAlias = Union[HighspyLeafTypes, Mapping[Any, 'HighspyNestedIndex'], Sequence[HighspyLeafTypes], np.ndarray[Any, np.dtype[Any]]]
    HighspyNestedResult: TypeAlias = Union[float, bool, Mapping[Any, 'HighspyNestedResult'], np.ndarray[Any, np.dtype[np.float64]]]
    HighspyIndexCollectionTypes: TypeAlias = Union[int, Integral, 'highs_var', 'highs_cons', 'highs_linear_expression', Mapping[Any, HighspyNestedIndex], Sequence[Union[int, Integral, 'highs_var', 'highs_cons', 'highs_linear_expression']], np.ndarray[Any, np.dtype[Any]]]
else:
    np_version = tuple((int(re.match('\\d+', part).group()) for part in np.__version__.split('.')))
    if sys.version_info >= (3, 9) and np_version >= (1, 22, 0):
        ndarray_object_type = np.ndarray[Any, np.dtype[np.object_]]
    else:
        ndarray_object_type = np.ndarray
    HighspyIndexCollectionTypes: TypeAlias = Union[int, Integral, 'highs_var', 'highs_cons', 'highs_linear_expression', Mapping[Any, Any], Sequence[Any], np.ndarray]
HighspyScalarTypes: TypeAlias = Union[float, int]
HighspyExpressionTypes: TypeAlias = 'highs_linear_expression'
HighspyLinearExpressionInputTypes: TypeAlias = Union[HighspyScalarTypes, 'highs_var', 'highs_linear_expression']
HighspyExpressionInputTypes: TypeAlias = Union['highs_var', 'highs_linear_expression']
HighspyConstraintTypes: TypeAlias = 'highs_linear_expression'
HighspyArrayItemTypes: TypeAlias = Union['highs_var', 'highs_linear_expression']

class Highs(_Highs):
    """
    HiGHS solver interface
    """
    __handle_keyboard_interrupt: bool = False
    __handle_user_interrupt: bool = False
    __solver_should_stop: bool = False
    __solver_stopped: RLock = RLock()
    __solver_started: Lock = Lock()
    __solver_status: Optional[HighsStatus] = None

    def __init__(self):
        super().__init__()
        self.callbacks = [HighsCallback(cb.HighsCallbackType(_), self) for _ in range(int(cb.HighsCallbackType.kCallbackMax) + 1)]
        self.enableCallbacks()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.is_solver_running():
            self.cancelSolve()
            self.wait()
        self.clear()

    def silent(self, turn_off_output: bool=True):
        """
        Disables solver output to the console.
        """
        super().setOptionValue('output_flag', not turn_off_output)

    def solve(self) -> Optional[HighsStatus]:
        """Runs the solver on the current problem.

        Returns:
            A HighsStatus object containing the solve status.
        """
        if not self.HandleKeyboardInterrupt:
            return super().run()
        else:
            return self.joinSolve(self.startSolve())

    def startSolve(self) -> Thread:
        """
        Starts the solver in a separate thread.  Useful for handling KeyboardInterrupts.
        Do not attempt to modify the model while the solver is running.

        Returns:
            A Thread object representing the solver thread.
        """
        if not self.is_solver_running():
            self.__solver_started.acquire()
            self.__solver_should_stop = False
            self.__solver_status = None
            t = Thread(target=Highs.__solve, args=(self,), daemon=True)
            t.start()
            try:
                self.__solver_started.acquire(True)
            finally:
                self.__solver_started.release()
            return t
        else:
            raise Exception('Solver is already running.')

    def is_solver_running(self) -> bool:
        is_running = True
        try:
            is_running = not self.__solver_stopped.acquire(False)
            return is_running
        finally:
            if not is_running:
                self.__solver_stopped.release()

    def __solve(self) -> None:
        try:
            self.__solver_stopped.acquire(True)
            self.__solver_started.release()
            self.__solver_status = super().run()
            _Highs.resetGlobalScheduler(False)
        finally:
            self.__solver_stopped.release()

    def joinSolve(self, solver_thread: Optional[Thread]=None, interrupt_limit: int=5) -> Optional[HighsStatus]:
        """
        Waits for the solver to finish. If solver_thread is provided, it will handle KeyboardInterrupts.

        Args:
            solver_thread: A Thread object representing the solver thread (optional).
            interrupt_limit: The number of times to allow KeyboardInterrupt before forcing termination (optional).

        Returns:
            A HighsStatus object containing the solve status.
        """
        result: Tuple[bool, Optional[HighsStatus]] = (False, None)
        if solver_thread is not None and interrupt_limit <= 0:
            try:
                while not result[0]:
                    result = self.wait(0.1)
                return result[1]
            except KeyboardInterrupt:
                print('KeyboardInterrupt: Waiting for HiGHS to finish...')
                self.cancelSolve()
        elif interrupt_limit > 0:
            for count in range(interrupt_limit):
                try:
                    while not result[0]:
                        result = self.wait(0.1)
                    return result[1]
                except KeyboardInterrupt:
                    print(f'Ctrl-C pressed {count + 1} times: Waiting for HiGHS to finish. ({interrupt_limit} times to force termination)')
                    self.cancelSolve()
            print('Forcing termination...')
            exit(1)
        try:
            self.__solver_stopped.acquire(True)
        except KeyboardInterrupt:
            pass
        finally:
            self.__solver_stopped.release()
        return self.__solver_status

    def wait(self, timeout: float=-1.0) -> Tuple[bool, Optional[HighsStatus]]:
        result: Tuple[bool, Optional[HighsStatus]] = (False, None)
        try:
            result = (self.__solver_stopped.acquire(True, timeout=timeout), self.__solver_status)
            return result
        finally:
            if result[0]:
                self.__solver_status = None
                self.__solver_stopped.release()

    def optimize(self) -> Optional[HighsStatus]:
        """
        Alias for the solve method.
        """
        return self.solve()

    def getObjective(self) -> Tuple[HighspyExpressionTypes, ObjSense]:
        """
        Retrieves the current objective function (as a linear expression) and sense.
        """
        lp = super().getLp()
        assert isinstance(lp.col_cost_, np.ndarray)
        nonzero_idxs = np.nonzero(lp.col_cost_)
        objective = highs_linear_expression()
        objective.idxs = nonzero_idxs[0].tolist()
        objective.vals = lp.col_cost_[nonzero_idxs].tolist()
        objective.constant = lp.offset_
        return (objective, super().getObjectiveSense()[1])

    def setObjective(self, obj: Optional[HighspyExpressionInputTypes]=None, sense: Optional[ObjSense]=None):
        """
        Updates the costs.

        Args:
            obj: An optional highs_linear_expression representing the new objective function.
            sense: An optional ObjSense value representing the new objective sense.

        Raises:
            Exception: If obj is an inequality or not a highs_linear_expression.
        """
        if obj is not None:
            expr = highs_linear_expression(obj) if isinstance(obj, highs_var) else obj
            if expr.bounds is not None:
                raise Exception('Objective cannot be an inequality')
            super().changeColsCost(self.numVariables, np.arange(self.numVariables, dtype=np.int32), np.full(self.numVariables, 0, dtype=np.float64))
            (idxs, vals) = expr.unique_elements()
            super().changeColsCost(len(idxs), idxs, vals)
            super().changeObjectiveOffset(expr.constant or 0.0)
        if sense is not None:
            super().changeObjectiveSense(sense)

    def minimize(self, obj: Optional[HighspyExpressionInputTypes]=None) -> Optional[HighsStatus]:
        """
        Solves a minimization of the objective and optionally updates the costs.

        Args:
            obj: An optional highs_linear_expression representing the new objective function.

        Raises:
            Exception: If obj is an inequality or not a highs_linear_expression.

        Returns:
            A HighsStatus object containing the solve status after minimization.
        """
        self.setObjective(obj, ObjSense.kMinimize)
        return self.solve()

    def maximize(self, obj: Optional[HighspyExpressionInputTypes]=None) -> Optional[HighsStatus]:
        """
        Solves a maximization of the objective and optionally updates the costs.

        Args:
            obj: An optional highs_linear_expression representing the new objective function.

        Raises:
            Exception: If obj is an inequality or not a highs_linear_expression.

        Returns:
            A HighsStatus object containing the solve status after maximization.
        """
        self.setObjective(obj, ObjSense.kMaximize)
        return self.solve()

    @staticmethod
    @overload
    def internal_get_value(array_values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]]], index_collection: Union[int, Integral, highs_var, highs_cons]) -> float:
        ...

    @staticmethod
    @overload
    def internal_get_value(array_values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]]], index_collection: highs_linear_expression) -> Union[float, bool]:
        ...

    @staticmethod
    @overload
    def internal_get_value(array_values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]]], index_collection: Mapping[Any, Mapping[Any, Union[int, Integral, highs_var, highs_cons]]]) -> Mapping[Any, Mapping[Any, float]]:
        ...

    @staticmethod
    @overload
    def internal_get_value(array_values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]]], index_collection: Mapping[Any, Union[int, Integral, highs_var, highs_cons]]) -> Mapping[Any, float]:
        ...

    @staticmethod
    @overload
    def internal_get_value(array_values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]]], index_collection: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @staticmethod
    @overload
    def internal_get_value(array_values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]]], index_collection: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    @staticmethod
    def internal_get_value(array_values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]]], index_collection: HighspyIndexCollectionTypes) -> Union[float, bool, Mapping[Any, Any], np.ndarray[Any, np.dtype[np.float64]]]:
        """
        Internal method to get the value of an index from an array of values. Could be value or dual, variable or constraint.
        """
        if isinstance(index_collection, (int, Integral, highs_var, highs_cons)):
            return array_values[int(index_collection)]
        elif isinstance(index_collection, highs_linear_expression):
            return index_collection.evaluate(array_values)
        elif isinstance(index_collection, Mapping):
            return {k: Highs.internal_get_value(array_values, v) for (k, v) in index_collection.items()}
        else:
            return np.asarray([Highs.internal_get_value(array_values, v) for v in index_collection])

    @overload
    def val(self, var: Union[int, Integral, highs_var, highs_cons]) -> float:
        ...

    @overload
    def val(self, var: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def val(self, var: Mapping[Any, Union[int, Integral, highs_var, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def val(self, var: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def val(self, var: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def val(self, var: HighspyIndexCollectionTypes):
        """
        Gets the value of a variable/index or expression in the solution.

        Args:
            var: A highs_var/index or highs_linear_expression object representing the variable.

        Returns:
            The value of the variable in the solution.
        """
        return Highs.internal_get_value(super().getSolution().col_value, var)

    @overload
    def vals(self, idxs: Union[int, Integral, highs_var, highs_cons]) -> float:
        ...

    @overload
    def vals(self, idxs: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def vals(self, idxs: Mapping[Any, Union[int, Integral, highs_var, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def vals(self, idxs: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def vals(self, idxs: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def vals(self, idxs: HighspyIndexCollectionTypes):
        """
        Gets the values of multiple variables in the solution.

        Args:
            idxs: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the solution values of the corresponding variables. If idxs is an iterable, returns a list of solution values for the variables.
        """
        return Highs.internal_get_value(super().getSolution().col_value, idxs)

    def variableName(self, var: Union[int, Integral, highs_var]):
        """
        Retrieves the name of a specific variable.

        Args:
            var: A highs_var object representing the variable.

        Raises:
            Exception: If the variable name cannot be found.

        Returns:
            The name of the specified variable.
        """
        [status, name] = super().getColName(int(var))
        failed = status != HighsStatus.kOk
        if failed:
            raise Exception('Variable name not found')
        return name

    @overload
    def variableNames(self, idxs: Mapping[Any, Union[highs_var, int, Integral]]) -> dict[Any, str]:
        ...

    @overload
    def variableNames(self, idxs: Iterable[Union[highs_var, int, Integral]]) -> list[str]:
        ...

    def variableNames(self, idxs: Iterable[Union[highs_var, int, Integral]]):
        """
        Retrieves the names of multiple variables.

        Args:
            idxs: An iterable of highs_var objects or a mapping where keys are identifiers and values are highs_var objects.

        Raises:
            Exception: If any variable name cannot be found.

        Returns:
            If idxs is a mapping, returns a dict where keys are the same keys from the input idxs and values are the names of the corresponding variables.
            If idxs is an iterable, returns a list of names for the specified variables.
        """
        if isinstance(idxs, Mapping):
            convert: Mapping[Any, Union[highs_var, int, Integral]] = idxs
            return {key: self.variableName(v) for (key, v) in convert.items()}
        else:
            return [self.variableName(v) for v in idxs]

    def allVariableNames(self) -> list[str]:
        """
        Retrieves the names of all variables in the model.

        Returns:
            A list of strings representing the names of all variables.
        """
        return super().getLp().col_names_

    @overload
    def variableValue(self, var: Union[int, Integral, highs_var]) -> float:
        ...

    @overload
    def variableValue(self, var: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def variableValue(self, var: Mapping[Any, Union[int, Integral, highs_var, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def variableValue(self, var: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def variableValue(self, var: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def variableValue(self, var: HighspyIndexCollectionTypes):
        """
        Retrieves the value of a specific variable in the solution.

        Args:
            var: A highs_var object representing the variable.

        Returns:
            The value of the specified variable in the solution.
        """
        return self.val(var)

    @overload
    def variableValues(self, idxs: Union[int, Integral, highs_var]) -> float:
        ...

    @overload
    def variableValues(self, idxs: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def variableValues(self, idxs: Mapping[Any, Union[int, Integral, highs_var, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def variableValues(self, idxs: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def variableValues(self, idxs: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def variableValues(self, idxs: HighspyIndexCollectionTypes):
        """
        Retrieves the values of multiple variables in the solution.

        Args:
            idxs: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the solution values of the corresponding variables. If idxs is an iterable, returns a list of solution values for the variables.
        """
        return self.vals(idxs)

    def allVariableValues(self):
        """
        Retrieves the values of all variables in the solution.

        Returns:
            A list of values for all variables in the solution.
        """
        return super().getSolution().col_value

    @overload
    def variableDual(self, var: Union[int, Integral, highs_var]) -> float:
        ...

    @overload
    def variableDual(self, var: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def variableDual(self, var: Mapping[Any, Union[int, Integral, highs_var, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def variableDual(self, var: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def variableDual(self, var: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def variableDual(self, var: HighspyIndexCollectionTypes):
        """
        Retrieves the dual value of a specific variable/index or expression in the solution.

        Args:
            var: A highs_var object representing the variable.

        Returns:
            The dual value of the specified variable in the solution.
        """
        return Highs.internal_get_value(super().getSolution().col_dual, var)

    @overload
    def variableDuals(self, idxs: Union[int, Integral, highs_var]) -> float:
        ...

    @overload
    def variableDuals(self, idxs: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def variableDuals(self, idxs: Mapping[Any, Union[int, Integral, highs_var, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def variableDuals(self, idxs: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def variableDuals(self, idxs: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def variableDuals(self, idxs: HighspyIndexCollectionTypes):
        """
        Retrieves the dual values of multiple variables in the solution.

        Args:
            idxs: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the dual values of the corresponding variables. If idxs is an iterable, returns a list of dual values for the variables.
        """
        return Highs.internal_get_value(super().getSolution().col_dual, idxs)

    def allVariableDuals(self):
        """
        Retrieves the dual values of all variables in the solution.

        Returns:
            A list of dual values for all variables in the solution.
        """
        return super().getSolution().col_dual

    @overload
    def constrValue(self, con: Union[int, Integral, highs_cons]) -> float:
        ...

    @overload
    def constrValue(self, con: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def constrValue(self, con: Mapping[Any, Union[int, Integral, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def constrValue(self, con: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def constrValue(self, con: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def constrValue(self, con: HighspyIndexCollectionTypes):
        """
        Retrieves the value of a specific constraint in the solution.

        Args:
            con: A highs_con object representing the constraint.

        Returns:
            The value of the specified constraint in the solution.
        """
        return Highs.internal_get_value(super().getSolution().row_value, con)

    @overload
    def constrValues(self, cons: Union[int, Integral, highs_cons]) -> float:
        ...

    @overload
    def constrValues(self, cons: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def constrValues(self, cons: Mapping[Any, Union[int, Integral, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def constrValues(self, cons: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def constrValues(self, cons: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def constrValues(self, cons: HighspyIndexCollectionTypes):
        """
        Retrieves the values of multiple constraints in the solution.

        Args:
            cons: A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.

        Returns:
            If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the solution values of the corresponding constraints. If cons is an iterable, returns a list of solution values for the constraints.
        """
        return Highs.internal_get_value(super().getSolution().row_value, cons)

    def allConstrValues(self):
        """
        Retrieves the values of all constraints in the solution.

        Returns:
            A list of values for all constraints in the solution.
        """
        return super().getSolution().row_value

    @overload
    def constrDual(self, con: Union[int, Integral, highs_cons]) -> float:
        ...

    @overload
    def constrDual(self, con: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def constrDual(self, con: Mapping[Any, Union[int, Integral, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def constrDual(self, con: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def constrDual(self, con: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def constrDual(self, con: HighspyIndexCollectionTypes):
        """
        Retrieves the dual value of a specific constraint in the solution.

        Args:
            con: A highs_con object representing the constraint.

        Returns:
            The dual value of the specified constraint in the solution.
        """
        return Highs.internal_get_value(super().getSolution().row_dual, con)

    @overload
    def constrDuals(self, cons: Union[int, Integral, highs_cons]) -> float:
        ...

    @overload
    def constrDuals(self, cons: highs_linear_expression) -> Union[float, bool]:
        ...

    @overload
    def constrDuals(self, cons: Mapping[Any, Union[int, Integral, highs_cons]]) -> Mapping[Any, float]:
        ...

    @overload
    def constrDuals(self, cons: Mapping[Any, Any]) -> Mapping[Any, Any]:
        ...

    @overload
    def constrDuals(self, cons: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]:
        ...

    def constrDuals(self, cons: HighspyIndexCollectionTypes):
        """
        Retrieves the dual values of multiple constraints in the solution.

        Args:
            cons: A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.

        Returns:
            If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the dual values of the corresponding constraints. If cons is an iterable, returns a list of dual values for the constraints.
        """
        return Highs.internal_get_value(super().getSolution().row_dual, cons)

    def allConstrDuals(self):
        """
        Retrieves the dual values of all constraints in the solution.

        Returns:
            A list of dual values for all constraints in the solution.
        """
        return super().getSolution().row_dual

    def addVariable(self, lb: float=0, ub: float=kHighsInf, obj: float=0.0, type: HighsVarType=HighsVarType.kContinuous, name: Optional[str]=None):
        """
        Adds a variable to the model.

        Args:
            lb: Lower bound of the variable (default is 0).
            ub: Upper bound of the variable (default is infinity).
            obj: Objective coefficient of the variable (default is 0).
            type: Type of the variable (continuous, integer; default is continuous).
            name: Optional name for the variable.

        Returns:
            A highs_var object representing the added variable.
        """
        status = super().addCol(obj, lb, ub, 0, np.empty(0, dtype=np.int32), np.empty(0, np.float64))
        if status != HighsStatus.kOk:
            raise Exception('Failed to add variable to the model.')
        var = highs_var(self.numVariables - 1, self)
        if type != HighsVarType.kContinuous:
            super().changeColIntegrality(var.index, type)
        if name is not None:
            super().passColName(var.index, name)
        return var

    @overload
    def addVariables(self, *nvars: int, out_array: Literal[True]=..., **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> HighspyArray:
        ...

    @overload
    def addVariables(self, *nvars: int, out_array: Literal[False], **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> dict[Any, highs_var]:
        ...

    @overload
    def addVariables(self, *nvars: Union[Mapping[Any, Any], Sequence[Any]], out_array: Literal[False]=..., **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> dict[Any, highs_var]:
        ...

    @overload
    def addVariables(self, *nvars: Union[Mapping[Any, Any], Sequence[Any]], out_array: Literal[True], **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> HighspyArray:
        ...

    @overload
    def addVariables(self, *nvars: Union[int, Mapping[Any, Any], Sequence[Any]], out_array: Optional[bool]=..., **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> Optional[Union[dict[Any, highs_var], HighspyArray]]:
        ...

    def addVariables(self, *nvars: Union[int, Mapping[Any, Any], Sequence[Any]], out_array: Optional[bool]=None, **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> Optional[Union[dict[Any, highs_var], HighspyArray]]:
        """
        Adds multiple variables to the model.

        Args:
            *args: A sequence of variables to be added. Can be a collection of scalars or indices (or mix).

            **kwargs: Optional keyword arguments.  Can be scalars, arrays, or mappings.
                lb: Lower bound of the variables (default is 0).
                ub: Upper bound of the variables (default is infinity).
                obj: Objective coefficient of the variables (default is 0).
                type: Type of the variables (continuous, integer; default is continuous).
                name: A collection of names for the variables (list or mapping).
                name_prefix: Prefix for the variable names.  Constructed name will be name_prefix + index.
                out_array: Return an array of highs_var objects instead of a dictionary.

        Returns:
            A highs_var collection (array or dictionary) representing the added variables.
        """
        if len(nvars) == 0:
            return None
        if out_array is not None and (not isinstance(out_array, bool)):
            raise TypeError(f'out_array expected bool or None, got {type(out_array).__name__}')
        shape = [n for n in nvars if isinstance(n, int)]
        expanded = [range(n) if isinstance(n, int) else n for n in nvars]
        indices: list[Any] = list(expanded[0] if len(expanded) == 1 else product(*expanded))
        N = len(indices)
        if len(shape) != len(nvars):
            shape = [N]

        def ensure_real(x: Union[Any, Mapping[Any, Any], Sequence[Any]]):
            if isinstance(x, (float, int)):
                return np.full(N, x, dtype=np.float64)
            elif isinstance(x, Mapping):
                mt: Mapping[Any, Any] = x
                if all((isinstance(v, (float, int)) for v in mt.values())):
                    m: Mapping[Any, Union[float, int]] = x
                    return np.fromiter((m[i] for i in indices), np.float64)
            elif isinstance(x, Sequence) and len(x) == N and all((isinstance(v, (float, int)) for v in x)):
                return np.asarray(x, dtype=np.float64)
            raise Exception('Invalid parameter.')

        def ensure_HighsVarType(x: Union[Any, Mapping[Any, Any], Sequence[Any]]):
            if x == HighsVarType.kContinuous:
                return None
            elif isinstance(x, HighsVarType):
                return np.full(N, x, dtype=np.uint8)
            elif isinstance(x, Mapping):
                mt: Mapping[Any, Any] = x
                if all((isinstance(v, HighsVarType) for v in mt.values())):
                    m: Mapping[Any, HighsVarType] = x
                    return np.fromiter((m[i] for i in indices), np.uint8)
            elif isinstance(x, Sequence) and len(x) == N and all((isinstance(v, HighsVarType) for v in x)):
                return np.asarray(x, dtype=np.uint8)
            raise Exception('Invalid parameter.')

        def ensure_optional_str(x: Optional[Any]):
            if x is None:
                return None
            elif isinstance(x, Sequence):
                mt: Sequence[Any] = x
                if len(mt) == N and all((isinstance(v, str) for v in mt)):
                    m: Sequence[str] = mt
                    return m
            raise Exception('Invalid parameter.')

        def ensure_str_or_none(x: Any) -> Optional[str]:
            if x is None or isinstance(x, str):
                return x
            else:
                raise Exception('Invalid parameter.')
        lb = ensure_real(kwargs.get('lb', 0.0))
        ub = ensure_real(kwargs.get('ub', kHighsInf))
        obj = ensure_real(kwargs.get('obj', 0))
        vartype = ensure_HighsVarType(kwargs.get('type', HighsVarType.kContinuous))
        name_prefix = ensure_str_or_none(kwargs.get('name_prefix', None))
        name = ensure_optional_str(kwargs.get('name', None))
        if out_array is None:
            out_array = all((isinstance(n, int) for n in nvars))
        start_idx = self.numVariables
        idx = np.arange(start_idx, start_idx + N, dtype=np.int32)
        status = super().addCols(N, obj, lb, ub, 0, np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32), np.empty(0, dtype=np.float64))
        if status != HighsStatus.kOk:
            raise Exception('Failed to add columns to the model.')
        if vartype is not None:
            super().changeColsIntegrality(N, idx, vartype)
        if name or name_prefix:
            names = name or [f'{name_prefix}{i}'.replace(' ', '') for i in indices]
            for (i, n) in zip(idx, names):
                super().passColName(int(i), str(n))
        return HighspyArray(np.asarray([highs_var(int(i), self) for i in idx]).reshape(shape), self) if out_array else {index: highs_var(int(i), self) for (index, i) in zip(indices, idx)}

    @overload
    def addIntegrals(self, *nvars: int, out_array: Optional[bool]=..., **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> HighspyArray:
        ...

    @overload
    def addIntegrals(self, *nvars: Union[Mapping[Any, Any], Sequence[Any]], out_array: Literal[False]=..., **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> dict[Any, highs_var]:
        ...

    @overload
    def addIntegrals(self, *nvars: Union[Mapping[Any, Any], Sequence[Any]], out_array: Literal[True], **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> HighspyArray:
        ...

    @overload
    def addIntegrals(self, *nvars: Union[int, Mapping[Any, Any], Sequence[Any]], out_array: Optional[bool]=..., **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> Optional[Union[dict[Any, highs_var], HighspyArray]]:
        ...

    def addIntegrals(self, *nvars: Union[int, Mapping[Any, Any], Sequence[Any]], out_array: Optional[bool]=None, **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]):
        """
        Alias for the addVariables method, for integer variables.
        """
        kwargs.setdefault('type', HighsVarType.kInteger)
        return self.addVariables(*nvars, out_array=out_array, **kwargs)

    @overload
    def addBinaries(self, *nvars: int, out_array: Optional[bool]=..., **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> HighspyArray:
        ...

    @overload
    def addBinaries(self, *nvars: Union[Mapping[Any, Any], Sequence[Any]], out_array: Literal[False]=..., **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> dict[Any, highs_var]:
        ...

    @overload
    def addBinaries(self, *nvars: Union[Mapping[Any, Any], Sequence[Any]], out_array: Literal[True], **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> HighspyArray:
        ...

    @overload
    def addBinaries(self, *nvars: Union[int, Mapping[Any, Any], Sequence[Any]], out_array: Optional[bool]=..., **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]) -> Optional[Union[dict[Any, highs_var], HighspyArray]]:
        ...

    def addBinaries(self, *nvars: Union[int, Mapping[Any, Any], Sequence[Any]], out_array: Optional[bool]=None, **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]]):
        """
        Alias for the addVariables method, for binary variables.
        """
        kwargs.setdefault('lb', 0)
        kwargs.setdefault('ub', 1)
        kwargs.setdefault('type', HighsVarType.kInteger)
        return self.addVariables(*nvars, out_array=out_array, **kwargs)

    def addIntegral(self, lb: float=0.0, ub: float=kHighsInf, obj: float=0.0, name: Optional[str]=None):
        """
        Alias for the addVariable method, for integer variables.
        """
        return self.addVariable(lb, ub, obj, HighsVarType.kInteger, name)

    def addBinary(self, obj: float=0.0, name: Optional[str]=None):
        """
        Alias for the addVariable method, for binary variables.
        """
        return self.addVariable(0, 1, obj, HighsVarType.kInteger, name)

    def deleteVariable(self, var_or_index: Union[int, Integral, highs_var, HighspyArrayItemTypes], *args: Union[Mapping[Any, highs_var], Iterable[Union[highs_var, HighspyArrayItemTypes]], highs_var, HighspyArray, HighspyArrayItemTypes]):
        """
        Deletes a variable from the model and updates the indices of subsequent variables in provided collections.

        Args:
            var_or_index: A highs_var object or an index representing the variable to be deleted.
            *args: Optional collections (lists, dicts, etc.) of highs_var objects whose indices need to be updated.
        """
        index = int(var_or_index)
        if index < self.numVariables:
            super().deleteVars(1, [index])
        for collection in args:
            if isinstance(collection, Mapping):
                mt: Mapping[Any, highs_var] = collection
                for (_, var) in mt.items():
                    if var.index > index:
                        var.index -= 1
                    elif var.index == index:
                        var.index = -1
            elif isinstance(collection, HighspyArray):
                for i in collection:
                    if isinstance(i, highs_var):
                        if i.index > index:
                            i.index -= 1
                        elif i.index == index:
                            i.index = -1
            elif isinstance(collection, Iterable):
                for var in collection:
                    if isinstance(var, highs_var):
                        if var.index > index:
                            var.index -= 1
                        elif var.index == index:
                            var.index = -1
            elif isinstance(collection, highs_var):
                if collection.index > index:
                    collection.index -= 1
                elif collection.index == index:
                    collection.index = -1

    def getVariables(self) -> list[highs_var]:
        """
        Retrieves all variables in the model.

        Returns:
            A list of highs_var objects, each representing a variable in the model.
        """
        return [highs_var(i, self) for i in range(self.numVariables)]

    @property
    def inf(self) -> float:
        """
        Represents infinity in the context of the solver.

        Returns:
            The value used to represent infinity.
        """
        return kHighsInf

    @property
    def numVariables(self) -> int:
        """
        Gets the number of variables in the model.

        Returns:
            The number of variables.
        """
        return super().getNumCol()

    @property
    def numConstrs(self) -> int:
        """
        Gets the number of constraints in the model.

        Returns:
            The number of constraints.
        """
        return super().getNumRow()

    def addConstr(self, expr: highs_linear_expression, name: Optional[str]=None) -> highs_cons:
        """
        Adds a constraint to the model.

        Args:
            expr: A highs_linear_expression to be added.
            name: Optional name of the constraint.

        Returns:
            A highs_cons object representing the added constraint.
        """
        con = self.__addRow(expr, self.numConstrs)
        if name is not None:
            super().passRowName(con.index, name)
        return con

    @overload
    def addConstrs(self, *args: highs_linear_expression, **kwargs: Optional[Union[str, Sequence[str]]]) -> list[highs_cons]:
        ...

    @overload
    def addConstrs(self, *args: Mapping[Any, highs_linear_expression], **kwargs: Optional[Union[str, Sequence[str]]]) -> dict[Any, highs_cons]:
        ...

    @overload
    def addConstrs(self, *args: Iterable[highs_linear_expression], **kwargs: Optional[Union[str, Sequence[str]]]) -> list[highs_cons]:
        ...

    @overload
    def addConstrs(self, *args: HighspyArray, **kwargs: Optional[Union[str, Sequence[str]]]) -> list[highs_cons]:
        ...

    def addConstrs(self, *args: Union[highs_linear_expression, Iterable[highs_linear_expression], HighspyArray], **kwargs: Optional[Union[str, Sequence[str]]]):
        """
        Adds multiple constraints to the model.

        Args:
            *args: A sequence of highs_linear_expression to be added.

            **kwargs: Optional keyword arguments.
                name_prefix: Prefix for the constraint names.  Constructed name will be name_prefix + index.
                name: A collection of names for the constraints (list or mapping).

        Returns:
            A highs_con collection array representing the added constraints.
        """
        name_prefix = kwargs.get('name_prefix', None)
        name = kwargs.get('name', None)
        generator = args[0] if len(args) == 1 and isinstance(args[0], Iterable) else args
        initial_rows = self.numConstrs
        cons: Union[dict[Any, highs_cons], list[highs_cons]]
        try:
            if isinstance(generator, Mapping):
                mt: Mapping[Any, highs_linear_expression] = generator
                cons = {key: self.__addRow(expr, initial_rows + count) for (count, (key, expr)) in enumerate(mt.items())}
            else:
                it: Iterable[highs_linear_expression] = cast(Iterable[highs_linear_expression], generator)
                cons = [self.__addRow(expr, initial_rows + count) for (count, expr) in enumerate(it)]
            if name or name_prefix:
                names = name or [f'{name_prefix}{n}' for n in range(self.numConstrs - initial_rows)]
                for (c, n) in zip(range(initial_rows, self.numConstrs), names):
                    super().passRowName(int(c), str(n))
        except Exception as e:
            status = super().deleteRows(self.numConstrs - initial_rows, np.arange(initial_rows, self.numConstrs, dtype=np.int32))
            if status != HighsStatus.kOk:
                raise Exception('Failed to rollback model after failure.  Model might be in a undeterminate state.')
            else:
                raise e
        return cons

    def __addRow(self, expr: highs_linear_expression, idx: int) -> highs_cons:
        """
        Internal method to add a constraint to the model.
        """
        if expr.bounds is not None:
            (idxs, vals) = expr.unique_elements()
            status = super().addRow(expr.bounds[0], expr.bounds[1], len(idxs), idxs, vals)
            if status != HighsStatus.kOk:
                raise Exception('Error adding constraint to the model.')
            return highs_cons(idx, self)
        else:
            raise Exception('Constraint bounds must be set via comparison (>=,==,<=).')

    def expr(self, optional: Optional[HighspyLinearExpressionInputTypes]=None) -> highs_linear_expression:
        """
        Creates a new highs_linear_expression object.

        Returns:
            A highs_linear_expression object.
        """
        return highs_linear_expression(optional)

    def getExpr(self, cons: Union[int, Integral, highs_cons]) -> highs_linear_expression:
        """
        Retrieves the highs_linear_expression of a constraint.

        Args:
            cons: A highs_con object or index representing the constraint.

        Returns:
            A highs_linear_expression object representing the expression of the constraint.
        """
        (status, lb, ub, nnz) = super().getRow(int(cons))
        if status != HighsStatus.kOk:
            raise Exception('Error retrieving constraint expression.')
        (status, idx, val) = super().getRowEntries(int(cons))
        if status != HighsStatus.kOk:
            raise Exception('Error retrieving constraint expression entries.')
        expr = highs_linear_expression()
        expr.bounds = (lb, ub)
        expr.idxs = list(idx)
        expr.vals = list(val)
        return expr

    def chgCoeff(self, cons: Union[highs_cons, int, Integral], var: Union[highs_var, int, Integral], val: float):
        """
        Changes the coefficient of a variable in a constraint.

        Args:
            cons: A highs_con object representing the constraint.
            var: A highs_var object representing the variable.
            val: The new coefficient value for the variable in the constraint.
        """
        super().changeCoeff(int(cons), int(var), val)

    def getConstrs(self) -> list[highs_cons]:
        """
        Retrieves all constraints in the model.

        Returns:
            A list of highs_cons objects, each representing a constraint in the model.
        """
        return [highs_cons(i, self) for i in range(self.numConstrs)]

    def removeConstr(self, cons_or_index: Union[highs_cons, int, Integral], *args: Union[Mapping[Any, highs_cons], Sequence[highs_cons], highs_cons]):
        """
        Removes a constraint from the model and updates the indices of subsequent constraints in provided collections.

        Args:
            cons_or_index: A highs_cons object or an index representing the constraint to be removed.
            *args: Optional collections (lists, dicts, etc.) of highs_cons objects whose indices need to be updated after the removal.
        """
        index = int(cons_or_index)
        if index < self.numConstrs:
            status = super().deleteRows(1, np.asarray([index], dtype=np.int32))
            if status != HighsStatus.kOk:
                raise Exception('Failed to delete constraint from the model.')
            for collection in args:
                if isinstance(collection, Mapping):
                    mt: Mapping[Any, highs_cons] = collection
                    for (_, con) in mt.items():
                        if con.index > index:
                            con.index -= 1
                        elif con.index == index:
                            con.index = -1
                elif isinstance(collection, Sequence):
                    for con in collection:
                        if con.index > index:
                            con.index -= 1
                        elif con.index == index:
                            con.index = -1
                elif collection.index > index:
                    collection.index -= 1
                elif collection.index == index:
                    collection.index = -1

    def setMinimize(self):
        """
        Sets the objective sense of the model to minimization.
        """
        super().changeObjectiveSense(ObjSense.kMinimize)

    def setMaximize(self):
        """
        Sets the objective sense of the model to maximization.
        """
        super().changeObjectiveSense(ObjSense.kMaximize)

    def setInteger(self, var_or_collection: Union[highs_var, int, HighspyArrayItemTypes, Iterable[Union[highs_var, int, HighspyArrayItemTypes]], HighspyArray]):
        """
        Sets a variable/collection to integer.

        Args:
            var_or_collection: A highs_var object/collection representing the variable to be set as integer.
        """
        if isinstance(var_or_collection, HighspyArray):
            idx = var_or_collection.idx()
            super().changeColsIntegrality(len(idx), idx, np.full(len(idx), HighsVarType.kInteger, dtype=np.uint8))
        elif isinstance(var_or_collection, Iterable):
            idx = np.fromiter(map(int, var_or_collection), dtype=np.int32)
            super().changeColsIntegrality(len(idx), idx, np.full(len(idx), HighsVarType.kInteger, dtype=np.uint8))
        else:
            super().changeColIntegrality(int(var_or_collection), HighsVarType.kInteger)

    def setContinuous(self, var_or_collection: Union[highs_var, int, HighspyArrayItemTypes, Iterable[Union[highs_var, int, HighspyArrayItemTypes]], HighspyArray]):
        """
        Sets a variable/collection to continuous.

        Args:
            var_or_collection: A highs_var object/collection representing the variable to be set as continuous.
        """
        if isinstance(var_or_collection, HighspyArray):
            idx = var_or_collection.idx()
            super().changeColsIntegrality(len(idx), idx, np.full(len(idx), HighsVarType.kContinuous, dtype=np.uint8))
        elif isinstance(var_or_collection, Iterable):
            idx = np.fromiter(map(int, var_or_collection), dtype=np.int32)
            super().changeColsIntegrality(len(idx), idx, np.full(len(idx), HighsVarType.kContinuous, dtype=np.uint8))
        else:
            super().changeColIntegrality(int(var_or_collection), HighsVarType.kContinuous)

    @staticmethod
    def idx(*args) -> np.ndarray[Any, np.dtype[np.int32]]:
        """Convert highs_var/highs_cons to a flat int32 index array.

        Can be called as:
            - ``h.idx(array)`` with a HighspyArray, numpy array, list, or tuple
            - ``h.idx(a, b, c)`` with individual highs_var or highs_cons objects

        Returns:
            A flat int32 numpy array of the underlying indices.
        """
        if len(args) == 1 and (not isinstance(args[0], (highs_var, highs_cons))):
            return np.asarray(args[0]).ravel().astype(np.int32)
        return np.array(args, dtype=np.int32).ravel()

    @staticmethod
    def qsum(items: Union[Iterable[HighspyArrayItemTypes], np.ndarray[Any, np.dtype[np.object_]]], initial: Optional[HighspyExpressionInputTypes]=None) -> HighspyExpressionTypes:
        """
        Performs a faster sum for highs_linear_expressions.

        Args:
            items: A collection of highs_linear_expressions or highs_vars to be summed.
        """
        expr = highs_linear_expression(initial)
        if isinstance(items, np.ndarray):
            X: np.ndarray[Any, np.dtype[np.object_]] = items
            for v in X.flat:
                expr += cast(Union[highs_var, highs_linear_expression], v)
        else:
            for item in items:
                expr += item
        return expr

    @staticmethod
    def __internal_callback(callback_type: cb.HighsCallbackType, message: str, data_out: cb.HighsCallbackOutput, data_in: Optional[cb.HighsCallbackInput], user_callback_data: Any):
        user_callback_data.callbacks[int(callback_type)].fire(callback_type, message, data_out, data_in)

    def enableCallbacks(self):
        """
        Enables callbacks, restarting them if they were previously enabled.
        """
        super().setCallback(Highs.__internal_callback, self)
        for c in self.callbacks:
            if len(c.callbacks) > 0:
                self.startCallback(c.callback_type)

    def clearCallbacks(self):
        """
        Clears all callbacks.
        """
        for c in self.callbacks:
            c.clear()

    def disableCallbacks(self):
        """
        Disables all callbacks, but does not clear them.
        """
        status = super().setCallback(None, None)
        if status != HighsStatus.kOk:
            raise Exception('Failed to disable callbacks.')

    def cancelSolve(self):
        """
        If HandleUserInterrupt is enabled, this method will signal the solver to stop.
        """
        self.__solver_should_stop = True

    @property
    def HandleKeyboardInterrupt(self) -> bool:
        """
        Get/Set whether the solver should handle KeyboardInterrupt (i.e., cancel solve on Ctrl+C). Also enables/disables HandleUserInterrupt.
        """
        return self.__handle_keyboard_interrupt

    @HandleKeyboardInterrupt.setter
    def HandleKeyboardInterrupt(self, value: bool):
        self.__handle_keyboard_interrupt = value
        self.HandleUserInterrupt = value

    @property
    def HandleUserInterrupt(self) -> bool:
        """
        Get/Set whether the solver should handle user interrupts (i.e., cancel solve on user request)
        """
        return self.__handle_user_interrupt

    @HandleUserInterrupt.setter
    def HandleUserInterrupt(self, value: bool):
        self.__handle_user_interrupt = value
        if value:
            self.cbSimplexInterrupt += self.__user_interrupt_event
            self.cbIpmInterrupt += self.__user_interrupt_event
            self.cbMipInterrupt += self.__user_interrupt_event
        else:
            self.cbSimplexInterrupt -= self.__user_interrupt_event
            self.cbIpmInterrupt -= self.__user_interrupt_event
            self.cbMipInterrupt -= self.__user_interrupt_event

    def __user_interrupt_event(self, e: HighsCallbackEvent):
        if self.__solver_should_stop:
            e.interrupt()

    class _callbackDescriptor:

        def __init__(self, callback_type: cb.HighsCallbackType):
            self.callback_type = callback_type

        def __get__(self, obj: Any, objtype: Any=None) -> Any:
            return self if obj is None else obj.callbacks[int(self.callback_type)]

        def __set__(self, obj: Any, value: Any) -> None:
            if obj.callbacks[int(self.callback_type)] is not value:
                raise Exception('Cannot set callback directly.  Use .subscribe(callback) instead.')
    cbLogging: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackLogging)
    cbSimplexInterrupt: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackSimplexInterrupt)
    cbIpmInterrupt: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackIpmInterrupt)
    cbMipSolution: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackMipSolution)
    cbMipImprovingSolution: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackMipImprovingSolution)
    cbMipLogging: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackMipLogging)
    cbMipInterrupt: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackMipInterrupt)
    cbMipGetCutPool: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackMipGetCutPool)
    cbMipDefineLazyConstraints: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackMipDefineLazyConstraints)
    cbMipUserSolution: HighsCallback = _callbackDescriptor(cb.HighsCallbackType.kCallbackMipUserSolution)

class HighsCallbackEvent(object):
    __slots__ = ['callback_type', 'message', 'data_out', 'data_in', 'user_data']

    def __init__(self, callback_type: cb.HighsCallbackType, message: str, data_out: cb.HighsCallbackOutput, data_in: Optional[cb.HighsCallbackInput], user_data: Optional[Any]):
        self.callback_type = callback_type
        self.message = message
        self.data_out = data_out
        self.data_in = data_in
        self.user_data = user_data

    def interrupt(self, interrupt_value: bool=True):
        """
        Sets the user interrupt flag in the callback data.
        """
        if self.data_in is not None:
            self.data_in.user_interrupt = interrupt_value

    def val(self, var_expr: HighspyIndexCollectionTypes):
        """
        Gets the value(s) of a variable/index or expression in the callback solution.
        """
        return Highs.internal_get_value(self.data_out.mip_solution, var_expr)

    def cut(self, index: int):
        """
        Gets the cut pool for the given index.
        """
        cut = highs_linear_expression()
        if self.data_out and index >= 0 and (index < self.data_out.cutpool_num_cut):
            (start, end) = self.data_out.cutpool_start[index:index + 2]
            cut.bounds = (self.data_out.cutpool_lower[index], self.data_out.cutpool_upper[index])
            cut.idxs = list(map(int, self.data_out.cutpool_index[start:end]))
            cut.vals = list(map(float, self.data_out.cutpool_value[start:end]))
        return cut

    @property
    def cuts(self):
        """
        Gets all cuts in the cut pool.
        """
        return [self.cut(i) for i in range(self.data_out.cutpool_num_cut)]

class HighsCallback(object):
    __slots__ = ['callbacks', 'user_callback_data', 'highs', 'callback_type']

    def __init__(self, callback_type: cb.HighsCallbackType, highs: Highs):
        self.callbacks: list[Callable[[HighsCallbackEvent], None]] = []
        self.user_callback_data: list[Any] = []
        self.callback_type = callback_type
        self.highs = proxy(highs)

    def subscribe(self, callback: Callable[[HighsCallbackEvent], None], user_data: Optional[Any]=None):
        """
        Subscribes a callback to the event.

        Args:
            callback: The callback function to be executed.
            user_data: Optional user data to be passed to the callback.
        """
        if len(self.callbacks) == 0:
            status = self.highs.startCallback(self.callback_type)
            if status != HighsStatus.kOk:
                raise Exception('Failed to start callback.')
        self.callbacks.append(callback)
        self.user_callback_data.append(user_data)
        return self

    def unsubscribe(self, callback: Callable[[HighsCallbackEvent], None]):
        """
        Unsubscribes a callback from the event.

        Args:
            callback: The callback function to be removed.
        """
        try:
            idx = self.callbacks.index(callback)
            del self.callbacks[idx]
            del self.user_callback_data[idx]
            if len(self.callbacks) == 0:
                self.highs.stopCallback(self.callback_type)
        except ValueError:
            pass
        return self

    def unsubscribe_by_data(self, user_data: Optional[Any]):
        """
        Unsubscribes a callback by user data.

        Args:
            user_data: The user data corresponding to the callback(s) to be removed.
        """
        idx = reversed([i for (i, ud) in enumerate(self.user_callback_data) if ud == user_data])
        for i in idx:
            del self.callbacks[i]
            del self.user_callback_data[i]
        if len(self.callbacks) == 0:
            self.highs.stopCallback(self.callback_type)
        return self

    def __iadd__(self, callback: Callable[[HighsCallbackEvent], None]):
        return self.subscribe(callback)

    def __isub__(self, callback: Callable[[HighsCallbackEvent], None]):
        return self.unsubscribe(callback)

    def clear(self):
        """
        Unsubscribes all callbacks from the event.
        """
        self.callbacks = []
        self.user_callback_data = []
        self.highs.stopCallback(self.callback_type)

    def fire(self, callback_type: cb.HighsCallbackType, message: str, data_out: cb.HighsCallbackOutput, data_in: cb.HighsCallbackInput):
        """
        Fires the event, executing all subscribed callbacks.
        """
        e = HighsCallbackEvent(callback_type, message, data_out, data_in, None)
        for (fn, user_data) in zip(self.callbacks, self.user_callback_data):
            e.user_data = user_data
            fn(e)

class HighspyArray(ndarray_object_type):
    """
    A numpy array wrapper for highs_var/highs_linear_expression objects.

    This provides additional type information for static analysis, and also allows faster sum operations.
    """
    highs: Optional[Highs]

    def __new__(cls, input_array: np.ndarray[Any, np.dtype[np.object_]], highs: Optional[Highs]) -> HighspyArray:
        obj = cast(HighspyArray, np.asarray(input_array).view(cls))
        obj.highs = highs
        return obj

    def __array_finalize__(self, obj: Optional[Any]):
        self.highs = getattr(obj, 'highs', None)

    @overload
    def __getitem__(self, key: Union[SupportsIndex, tuple[SupportsIndex, ...]]) -> HighspyArrayItemTypes:
        ...

    @overload
    def __getitem__(self, key: Union[slice, Sequence[int], np.ndarray[Any, np.dtype[np.integer[Any]]], np.ndarray[Any, np.dtype[np.bool_]], tuple[Union[None, slice, SupportsIndex, np.ndarray[Any, np.dtype[Any]]], ...]]) -> HighspyArray:
        ...

    @overload
    def __getitem__(self, key: Any) -> Union[HighspyArray, HighspyArrayItemTypes]:
        ...

    def __getitem__(self, key: Any) -> Union[HighspyArray, HighspyArrayItemTypes]:
        return super(HighspyArray, self).__getitem__(key)

    def __iter__(self) -> Iterator[HighspyArrayItemTypes]:
        return super().__iter__()

    def __ge__(self, other: Any) -> HighspyArray:
        return cast(HighspyArray, np.greater_equal(self, other, dtype=np.object_))

    def __le__(self, other: Any) -> HighspyArray:
        return cast(HighspyArray, np.less_equal(self, other, dtype=np.object_))

    def __eq__(self, other: Any) -> HighspyArray:
        return cast(HighspyArray, np.equal(self, other, dtype=np.object_))

    @overload
    def sum(self, axis: None=None, dtype: Optional[Any]=None, out: None=...) -> highs_linear_expression:
        ...

    @overload
    def sum(self, axis: Any, dtype: Optional[Any]=None, out: HighspyArray=...) -> HighspyArray:
        ...

    def sum(self, axis: Optional[int]=None, dtype: Optional[Any]=None, out: Optional[np.ndarray[Any, np.dtype[np.object_]]]=None, **unused_kwargs: Any) -> Union[HighspyArray, highs_linear_expression]:
        if self.highs is not None:
            if axis is not None:
                return HighspyArray(np.apply_along_axis(self.highs.qsum, axis, self, initial=unused_kwargs.get('initial', None)), self.highs)
            else:
                return self.highs.qsum(self, unused_kwargs.get('initial', None))
        else:
            raise Exception('Cannot sum without a Highs object.')

    def idx(self) -> np.ndarray[Any, np.dtype[np.int32]]:
        """Convert to a flat int32 index array for passing to the HiGHS C++ API.

        Each element's ``__index__`` method is called to extract its integer
        index (e.g., ``highs_var.index`` or ``highs_cons.index``).
        The result is always 1-D, regardless of the array's shape.

        Returns:
            A new flat int32 numpy array of the underlying indices.
        """
        return self.ravel().astype(np.int32)

class highs_var(object):
    """
    Variable index wrapper for HiGHS
    """
    __slots__ = ['index', 'highs']

    def __init__(self, i: int, highs: Highs):
        self.index = i
        self.highs = proxy(highs)

    def __repr__(self):
        return f'highs_var({self.index})'

    @property
    def name(self) -> str:
        return self.highs.variableName(self)

    @name.setter
    def name(self, value: str):
        self.highs.passColName(self.index, value)

    def __int__(self):
        return int(self.index)

    def __index__(self) -> int:
        return int(self.index)

    def __hash__(self):
        return int(self.index)

    def __neg__(self):
        expr = highs_linear_expression()
        expr.idxs = [self.index]
        expr.vals = [-1.0]
        return expr

    def __add__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr += other
        return expr

    def __radd__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr += other
        return expr

    def __mul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr *= other
        return expr

    def __rmul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr *= other
        return expr

    def __rsub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(other)
        expr -= self
        return expr

    def __sub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr -= other
        return expr

    def __le__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, highs_linear_expression):
            return other.__ge__(self)
        else:
            return highs_linear_expression(self).__le__(other)

    @overload
    def __eq__(self, other: None) -> bool:
        ...

    @overload
    def __eq__(self, other: HighspyLinearExpressionInputTypes) -> highs_linear_expression:
        ...

    def __eq__(self, other: Any) -> Union[bool, highs_linear_expression]:
        if other is None:
            return True
        elif isinstance(other, highs_linear_expression):
            return other.__eq__(self)
        else:
            return highs_linear_expression(self).__eq__(other)

    def __ne__(self, other: Optional[Any]) -> bool:
        if other is None:
            return True
        else:
            raise Exception('Invalid comparison.')

    def __ge__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, highs_linear_expression):
            return other.__le__(self)
        else:
            return highs_linear_expression(self).__ge__(other)

    def __truediv__(self, other: Union[float, int, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr /= other
        return expr

    def __rtruediv__(self, other: Union[float, int, highs_linear_expression]):
        raise Exception('Only division of a linear expression by a scalar is allowed.')

class highs_cons(object):
    """
    Constraint index wrapper for HiGHS
    """
    __slots__ = ['index', 'highs']

    def __init__(self, i: int, highs: Highs):
        self.index = i
        self.highs = proxy(highs)

    def __repr__(self):
        return f'highs_cons({self.index})'

    def __int__(self):
        return int(self.index)

    def __index__(self) -> int:
        return int(self.index)

    def __hash__(self):
        return int(self.index)

    def expr(self) -> highs_linear_expression:
        """
        Retrieves the expression of the constraint.
        """
        return self.highs.getExpr(self)

    @property
    def name(self) -> str:
        (status, name) = self.highs.getRowName(self.index)
        if status != HighsStatus.kOk:
            raise Exception('Error retrieving constraint name.')
        return name

    @name.setter
    def name(self, value: str):
        self.highs.passRowName(self.index, value)

class highs_linear_expression(object):
    """
    Linear constraint builder for HiGHS
    """
    __slots__ = ['idxs', 'vals', 'constant', 'bounds']

    @overload
    def __init__(self, other: None=None) -> None:
        ...

    @overload
    def __init__(self, other: Union[float, int]) -> None:
        ...

    @overload
    def __init__(self, other: highs_var) -> None:
        ...

    @overload
    def __init__(self, other: highs_linear_expression) -> None:
        ...

    def __init__(self, other: Optional[Union[float, int, highs_var, highs_linear_expression]]=None):
        self.constant: Optional[float] = None
        self.bounds: Optional[tuple[float, float]] = None
        if other is None:
            self.idxs: list[int] = []
            self.vals: list[float] = []
        elif isinstance(other, highs_linear_expression):
            self.idxs = list(other.idxs)
            self.vals = list(other.vals)
            self.constant = other.constant
            self.bounds = other.bounds if other.bounds is not None else None
        elif isinstance(other, highs_var):
            self.idxs = [other.index]
            self.vals = [1.0]
        else:
            self.idxs = []
            self.vals = []
            self.constant = float(other)

    def simplify(self):
        """
        Simplifies the linear expression by combining duplicate variables.
        """
        copy = highs_linear_expression()
        (copy.idxs, copy.vals) = (v.tolist() for v in self.unique_elements())
        copy.bounds = self.bounds if self.bounds is not None else None
        copy.constant = self.constant
        return copy

    def copy(self):
        """
        Creates a copy of the linear expression.
        """
        return highs_linear_expression(self)

    def evaluate(self, values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]]]) -> Union[float, bool]:
        """
        Evaluates the linear expression given a solution array (values).
        """
        result = sum((v * values[c] for (c, v) in zip(self.idxs, self.vals))) + (self.constant or 0.0)
        return result if self.bounds is None else self.bounds[0] <= result <= self.bounds[1]

    def __repr__(self):
        v = str.join('  ', [f'{c}_v{x}' for (x, c) in zip(self.idxs, self.vals)])
        if self.bounds is None:
            return f'{v}' + (f'  {self.constant}' if self.constant is not None else '')
        elif self.bounds[0] == self.bounds[1]:
            return f'{v} == {self.bounds[0] - (self.constant or 0.0)}'
        else:
            return f'{self.bounds[0]} <= {v} <= {self.bounds[1]}'

    def __str__(self):
        (idxs, vals) = self.unique_elements()
        v = str.join('  ', [f'{c}_v{x}' for (x, c) in zip(idxs, vals)])
        if self.bounds is None:
            return f'{v}' + (f'  {self.constant}' if self.constant is not None else '')
        elif self.bounds[0] == self.bounds[1]:
            return f'{v} == {self.bounds[0] - (self.constant or 0.0)}'
        else:
            return f'{self.bounds[0]} <= {v} <= {self.bounds[1]}'

    def __ne__(self, other: Optional[Any]) -> bool:
        if other is None:
            return True
        else:
            raise Exception('Invalid comparison.')

    @overload
    def __eq__(self, other: None) -> bool:
        ...

    @overload
    def __eq__(self, other: HighspyLinearExpressionInputTypes) -> highs_linear_expression:
        ...

    @overload
    def __eq__(self, other: Sequence[Union[float, int]]) -> highs_linear_expression:
        ...

    def __eq__(self, other: Any) -> Union[bool, highs_linear_expression]:
        if self.bounds is not None:
            raise Exception('Bounds have already been set.')
        elif isinstance(other, (float, int)):
            copy = highs_linear_expression(self)
            copy.bounds = (float(other) - (self.constant or 0.0), float(other) - (self.constant or 0.0))
            copy.constant = None
            return copy
        elif isinstance(other, highs_linear_expression):
            if other.bounds is not None:
                raise Exception('Bounds have already been set.')
            copy = highs_linear_expression()
            if len(other.idxs) > len(self.idxs) or (len(other.idxs) == len(self.idxs) and other.constant is None and (self.constant is not None)):
                copy.idxs = other.idxs + self.idxs
                copy.vals = other.vals + [-v for v in self.vals]
                copy.bounds = ((self.constant or 0.0) - (other.constant or 0.0), (self.constant or 0.0) - (other.constant or 0.0))
            else:
                copy.idxs = self.idxs + other.idxs
                copy.vals = self.vals + [-v for v in other.vals]
                copy.bounds = ((other.constant or 0.0) - (self.constant or 0.0), (other.constant or 0.0) - (self.constant or 0.0))
            return copy
        elif isinstance(other, highs_var):
            copy = highs_linear_expression()
            if len(self.idxs) == 0 or (len(self.idxs) == 1 and self.constant is not None):
                copy.idxs = [other.index] + self.idxs
                copy.vals = [1.0] + [-v for v in self.vals]
                copy.bounds = (self.constant or 0.0, self.constant or 0.0)
            else:
                copy.idxs = self.idxs + [other.index]
                copy.vals = self.vals + [-1.0]
                copy.bounds = (-(self.constant or 0.0), -(self.constant or 0.0))
            return copy
        elif hasattr(other, '__getitem__') and hasattr(other, '__len__') and (len(other) == 2):
            if not (isinstance(other[0], (float, int)) and isinstance(other[1], (float, int))):
                raise Exception('Provided bounds were not valid numbers.')
            copy = highs_linear_expression(self)
            copy.bounds = (float(other[0]) - (copy.constant or 0.0), float(other[1]) - (copy.constant or 0.0))
            copy.constant = None
            return copy
        else:
            raise Exception('Unknown comparison.')

    def __le__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if self.bounds is not None:
            raise Exception('Bounds have already been set.')
        elif self.__is_active_chain():
            other = other if isinstance(other, highs_linear_expression) else highs_linear_expression(other)
            order = self.__get_chain(other, False)
            if self.__is_equal_except_bounds(order[2]):
                return self.__compose_chain(*order[1:])
            elif other.__is_equal_except_bounds(order[1]):
                return self.__compose_chain(*order[:-1])
        self.__reset_chain(self, other, None)
        if isinstance(other, (float, int)):
            copy = highs_linear_expression(self)
            copy.bounds = (-kHighsInf, float(other) - (copy.constant or 0.0))
            copy.constant = None
            return copy
        elif isinstance(other, highs_linear_expression):
            if other.bounds is None:
                copy = highs_linear_expression()
                if len(other.idxs) > len(self.idxs) or (len(other.idxs) == len(self.idxs) and other.constant is None and (self.constant is not None)):
                    copy.idxs = other.idxs + self.idxs
                    copy.vals = other.vals + [-v for v in self.vals]
                    copy.bounds = ((self.constant or 0.0) - (other.constant or 0.0), kHighsInf)
                else:
                    copy.idxs = self.idxs + other.idxs
                    copy.vals = self.vals + [-v for v in other.vals]
                    copy.bounds = (-kHighsInf, (other.constant or 0.0) - (self.constant or 0.0))
                return copy
            else:
                raise Exception('Bounds have already been set.')
        else:
            copy = highs_linear_expression()
            if len(self.idxs) == 0 or (len(self.idxs) == 1 and self.constant is not None):
                copy.idxs = [other.index] + self.idxs
                copy.vals = [1.0] + [-v for v in self.vals]
                copy.bounds = (self.constant or 0.0, kHighsInf)
            else:
                copy.idxs = self.idxs + [other.index]
                copy.vals = self.vals + [-1.0]
                copy.bounds = (-kHighsInf, -(self.constant or 0.0))
            return copy

    def __ge__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if self.bounds is not None:
            raise Exception('Bounds have already been set.')
        elif self.__is_active_chain():
            other = other if isinstance(other, highs_linear_expression) else highs_linear_expression(other)
            order = self.__get_chain(other, True)
            if self.__is_equal_except_bounds(order[1]):
                return self.__compose_chain(*order[:-1])
            elif other.__is_equal_except_bounds(order[2]):
                return self.__compose_chain(*order[1:])
        self.__reset_chain(None, other, self)
        if isinstance(other, (float, int)):
            copy = highs_linear_expression(self)
            copy.bounds = (float(other) - (self.constant or 0.0), kHighsInf)
            copy.constant = None
            return copy
        elif isinstance(other, highs_linear_expression):
            if other.bounds is None:
                copy = highs_linear_expression()
                if len(self.idxs) > len(other.idxs) or (len(self.idxs) == len(other.idxs) and self.constant is None and (other.constant is not None)):
                    copy.idxs = self.idxs + other.idxs
                    copy.vals = self.vals + [-v for v in other.vals]
                    copy.bounds = ((other.constant or 0.0) - (self.constant or 0.0), kHighsInf)
                else:
                    copy.idxs = other.idxs + self.idxs
                    copy.vals = other.vals + [-v for v in self.vals]
                    copy.bounds = (-kHighsInf, (self.constant or 0.0) - (other.constant or 0.0))
                return copy
            else:
                raise Exception('Bounds have already been set.')
        else:
            copy = highs_linear_expression()
            if len(self.idxs) > 1:
                copy.idxs = self.idxs + [other.index]
                copy.vals = self.vals + [-1.0]
                copy.bounds = (-(self.constant or 0.0), kHighsInf)
            else:
                copy.idxs = [other.index] + self.idxs
                copy.vals = [1.0] + [-v for v in self.vals]
                copy.bounds = (-kHighsInf, self.constant or 0.0)
            return copy

    def __radd__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        return self + other

    def __add__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        copy = highs_linear_expression(self)
        copy += other
        return copy

    def __neg__(self):
        copy = highs_linear_expression()
        copy.idxs = list(self.idxs)
        copy.vals = [-v for v in self.vals]
        copy.constant = -self.constant if self.constant is not None else None
        if self.bounds is not None:
            copy.bounds = (-self.bounds[1], -self.bounds[0])
        return copy

    def __rmul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        return self * other

    def __mul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        copy = highs_linear_expression(self)
        copy *= other
        return copy

    def __rsub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        copy = highs_linear_expression(other)
        copy -= self
        return copy

    def __sub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        copy = highs_linear_expression(self)
        copy -= other
        return copy

    def __int__(self):
        raise TypeError('Cannot convert to int')

    def unique_elements(self):
        """
        Collects unique variables and sums their corresponding values.  Keeps all values (including zeros).
        """
        groups = np.asarray(self.idxs, dtype=np.int32)
        order = np.argsort(groups, kind='stable')
        groups = groups[order]
        index = np.ones(len(groups), dtype=bool)
        index[:-1] = groups[1:] != groups[:-1]
        if index.all():
            values = np.asarray(self.vals, dtype=np.float64)
            values = values[order]
            return (groups, values)
        else:
            values = np.asarray(self.vals, dtype=np.float64)
            values = np.cumsum(values[order])
            values = values[index]
            groups = groups[index]
            values[1:] = values[1:] - values[:-1]
            return (groups, values)

    def reduced_elements(self):
        """
        Similar to unique_elements, except keeps only non-zero values
        """
        (vx, vl) = self.unique_elements()
        zx = np.nonzero(vl)
        return (vx[zx], vl[zx])

    def __iadd__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, highs_var):
            self.idxs.append(other.index)
            self.vals.append(1.0)
            return self
        elif isinstance(other, highs_linear_expression):
            if self.constant is not None and other.bounds is not None or (self.bounds is not None and other.constant is not None):
                raise Exception('Cannot add a bounded constraint to a constraint with a constant, i.e., (lb <= expr1 <= ub) + (expr2 + c). \n                    Unsure of your intent. Did you want: lb + c <= expr1 + expr2 <= ub + c?  Try: (lb <= expr1 <= ub) + (expr2 == c) instead.')
            self.idxs.extend(other.idxs)
            self.vals.extend(other.vals)
            if self.constant is not None or other.constant is not None:
                self.constant = (self.constant or 0.0) + (other.constant or 0.0)
            if self.bounds is not None and other.bounds is not None:
                self.bounds = (self.bounds[0] + other.bounds[0], self.bounds[1] + other.bounds[1])
            elif self.bounds is None and other.bounds is not None:
                self.bounds = other.bounds
            return self
        else:
            if self.bounds is not None:
                raise Exception('Cannot add a constant to a bounded constraint, i.e., (lb <= expr <= ub) + c. \n                    Unsure of your intent. Did you want: lb + c <= expr <= ub + c?  Try: (lb <= expr <= ub) + (highs_linear_expression() == c) instead.')
            self.constant = float(other) + (self.constant or 0.0)
            return self

    def __isub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, highs_var):
            self.idxs.append(other.index)
            self.vals.append(-1.0)
            return self
        elif isinstance(other, highs_linear_expression):
            if self.constant is not None and other.bounds is not None or (self.bounds is not None and other.constant is not None):
                raise Exception('Cannot subtract a bounded constraint to a constraint with a constant, i.e., (lb <= expr1 <= ub) - (expr2 + c). \n                    Unsure of your intent. Did you want: lb - c <= expr1 - expr2 <= ub - c?  Try: (lb <= expr1 <= ub) - (expr2 == c) instead.')
            self.idxs.extend(other.idxs)
            self.vals.extend([-v for v in other.vals])
            if self.constant is not None or other.constant is not None:
                self.constant = (self.constant or 0.0) - (other.constant or 0.0)
            if self.bounds is not None and other.bounds is not None:
                self.bounds = (self.bounds[0] - other.bounds[1], self.bounds[1] - other.bounds[0])
            elif self.bounds is None and other.bounds is not None:
                self.bounds = (-other.bounds[1], -other.bounds[0])
            return self
        else:
            if self.bounds is not None:
                raise Exception('Cannot subtract a constant to a bounded constraint, i.e., (lb <= expr <= ub) - c. \n                    Unsure of your intent. Did you want: lb - c <= expr <= ub - c?  Try: (lb <= expr <= ub) - (highs_linear_expression() == c) instead.')
            self.constant = (self.constant or 0.0) - float(other)
            return self

    def __truediv__(self, other: Union[float, int, highs_linear_expression]):
        copy = highs_linear_expression(self)
        copy /= other
        return copy

    def __itruediv__(self, other: Union[float, int, highs_linear_expression]):
        if isinstance(other, (float, int)):
            divisor = float(other)
        elif isinstance(other, highs_linear_expression) and other.idxs == [] and (other.constant is not None):
            divisor = float(other.constant)
        else:
            raise Exception('Only division by a scalar is allowed.')
        if divisor == 0:
            raise ZeroDivisionError('division by zero')
        return self.__imul__(1.0 / divisor)

    def __rtruediv__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        raise Exception('Only division of a linear expression by a scalar is allowed.')

    def __imul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, (float, int)):
            scale = float(other)
        elif isinstance(other, highs_linear_expression) and other.idxs == [] and (other.constant is not None):
            scale = float(other.constant)
        else:
            scale = None
        if scale is not None:
            self.vals = [scale * v for v in self.vals]
            if self.constant is not None:
                self.constant *= scale
            if self.bounds is not None:
                if scale >= 0:
                    self.bounds = (scale * self.bounds[0], scale * self.bounds[1])
                else:
                    self.bounds = (scale * self.bounds[1], scale * self.bounds[0])
            return self
        elif self.idxs == [] and self.constant is not None:
            scale = self.constant
            if isinstance(other, highs_linear_expression):
                self.idxs = other.idxs
                self.vals = [scale * v for v in other.vals]
                if other.constant is not None:
                    self.constant = other.constant * scale
                else:
                    self.constant = None
                if other.bounds is not None:
                    if scale >= 0:
                        self.bounds = (scale * other.bounds[0], scale * other.bounds[1])
                    else:
                        self.bounds = (scale * other.bounds[1], scale * other.bounds[0])
            elif isinstance(other, highs_var):
                self.idxs = [other.index]
                self.vals = [scale]
                self.constant = None
            else:
                raise Exception('Unexpected parameters.')
            return self
        elif isinstance(other, highs_var):
            raise Exception('Only linear expressions are allowed.')
        else:
            raise Exception('Unexpected parameters.')
    __chain = local()

    def __bool__(self):
        highs_linear_expression.__chain.check = self if self.bounds is not None and (self.bounds[0] != self.bounds[1] or (self.bounds[0] == kHighsInf or self.bounds[1] == -kHighsInf)) else None
        LHS = getattr(highs_linear_expression.__chain, 'left', None)
        EXR = getattr(highs_linear_expression.__chain, 'inner', None)
        RHS = getattr(highs_linear_expression.__chain, 'right', None)
        highs_linear_expression.__chain.left = highs_linear_expression(LHS) if LHS is not None else None
        highs_linear_expression.__chain.inner = highs_linear_expression(EXR) if EXR is not None else None
        highs_linear_expression.__chain.right = highs_linear_expression(RHS) if RHS is not None else None
        return True

    def __is_equal_except_bounds(self, other: highs_linear_expression) -> bool:
        return self.idxs == other.idxs and self.vals == other.vals and (self.constant == other.constant)

    def __is_active_chain(self):
        return getattr(highs_linear_expression.__chain, 'check', None) is not None

    def __reset_chain(self, LHS: Optional[Any]=None, inner: Optional[Any]=None, RHS: Optional[Any]=None) -> None:
        highs_linear_expression.__chain.check = None
        highs_linear_expression.__chain.left = LHS
        highs_linear_expression.__chain.inner = inner
        highs_linear_expression.__chain.right = RHS

    def __get_chain(self, other: highs_linear_expression, is_ge_than: bool):
        LHS = getattr(highs_linear_expression.__chain, 'left', None)
        RHS = getattr(highs_linear_expression.__chain, 'right', None)
        inner = getattr(highs_linear_expression.__chain, 'inner', None)
        assert (LHS is None) ^ (RHS is None) == 1
        order = np.asarray([expr for expr in [other, LHS, inner, RHS, self] if expr is not None])
        if not is_ge_than:
            (order[0], order[-1]) = (order[-1], order[0])
        return order

    def __compose_chain(self, left: highs_linear_expression, inner: highs_linear_expression, right: highs_linear_expression):
        self.__reset_chain()
        assert not isinstance(inner, highs_linear_expression) or inner.bounds is None, 'Bounds already set in chain comparison.'
        (LHS_vars, LHS_vals) = left.reduced_elements()
        (RHS_vars, RHS_vals) = right.reduced_elements()
        if not np.array_equal(LHS_vars, RHS_vars) or not np.array_equal(LHS_vals, RHS_vals):
            raise Exception('Mismatched variables in chain comparison.')
        if len(LHS_vars) > len(inner.idxs):
            copy = highs_linear_expression()
            copy.idxs = LHS_vars.tolist() + inner.idxs
            copy.vals = LHS_vals.tolist() + [-v for v in inner.vals]
            copy.bounds = ((inner.constant or 0.0) - (right.constant or 0.0), (inner.constant or 0.0) - (left.constant or 0.0))
        else:
            copy = highs_linear_expression(inner)
            copy.idxs.extend(LHS_vars)
            copy.vals.extend([-v for v in LHS_vals])
            copy.bounds = ((left.constant or 0.0) - (copy.constant or 0.0), (right.constant or 0.0) - (copy.constant or 0.0))
            copy.constant = None
        return copy

def qsum(items: Iterable[HighspyArrayItemTypes], initial: Optional[HighspyExpressionInputTypes]=None):
    """
    Performs a faster sum for highs_linear_expressions.

    Args:
        items: A collection of highs_linear_expressions or highs_vars to be summed.
    """
    return Highs.qsum(items, initial)
